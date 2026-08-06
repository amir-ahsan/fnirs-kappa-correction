#!/usr/bin/env python3
"""
LEGACY / NOT USED FOR THE MANUSCRIPT TABLES.

This early isotropic (g=0) CSF Monte Carlo has been SUPERSEDED by the anisotropic
(Henyey-Greenstein, g=0.9) transport in mc_2layer.py / mc_production.py, which is
what produces every reported f_cortex and CSF ratio. mc_csf.py is retained only
for historical reference; it is not imported by any analysis script.

Vectorized white-Monte-Carlo for the layered fNIRS slab, to compute the
cortical sensitivity fraction f_cortex = <L_cortex> / <L_total>, where
<L_layer> is the mean partial pathlength in a layer over detected photons
(weighted by Beer-Lambert survival w = exp(-sum mu_a*L)).  This is the
Monte-Carlo analogue of the adjoint/Jacobian sensitivity fraction
S_cortex/S_total used in the paper.

Scattering: isotropic (g=0) with mu_s = mu_s' (the reduced-scattering
"similarity" MC, appropriate when only mu_s' is specified).
Surface: index-mismatched (n_tissue=1.4 / n_air=1.0) with Fresnel reflectance.
Detection: photons exiting the top surface are binned by exit radius, so
every exiting photon contributes to some SDS annulus (efficient).

Outputs per SDS bin: sum(w*L_cortex), sum(w*L_total), sum(w), count.
Accumulators are additive across independent runs (different seeds).
"""
import numpy as np, json, time

def run(geometry, N, seed, sds_centers, half_width=2.5,
        z_max=70.0, max_iter=4000, n_tissue=1.4):
    rng = np.random.default_rng(seed)
    bounds = np.array(geometry['bounds'], float)        # layer lower-z boundaries, e.g. [0,10,12]
    mua    = np.array(geometry['mua'], float)
    mus    = np.array(geometry['musp'], float)          # g=0 => mu_s = mu_s'
    nlay   = len(mua)
    cortex_z = geometry['cortex_z']                     # z above which = cortex
    bounds_ext = np.append(bounds, z_max)               # lower edge of each layer incl. cap

    # photon state
    x = np.zeros(N); y = np.zeros(N); z = np.zeros(N)
    ux = np.zeros(N); uy = np.zeros(N); uz = np.ones(N) # launch straight down (+z)
    L = np.zeros((N, nlay))                             # partial pathlength per layer
    layer = np.zeros(N, int)
    active = np.ones(N, bool)
    detected = np.zeros(N, bool)
    exit_r = np.full(N, -1.0)


    for it in range(max_iter):
        idx = np.where(active)[0]
        if idx.size == 0:
            break
        lay = layer[idx]
        mus_i = mus[lay]; zz = z[idx]; uzz = uz[idx]
        # step length (white MC: scattering only)
        step = -np.log(rng.random(idx.size)) / mus_i
        # distance to the z-boundary in the direction of travel
        lower = bounds_ext[lay + 1]      # deeper boundary of current layer (or z_max cap)
        upper = bounds[lay]              # shallower boundary (top); for layer 0 this is 0 (surface)
        target = np.where(uzz > 0, lower, upper)
        with np.errstate(divide='ignore', invalid='ignore'):
            dist_b = (target - zz) / uzz
        dist_b = np.where(np.abs(uzz) < 1e-12, np.inf, dist_b)
        dist_b = np.where(dist_b < 0, np.inf, dist_b)    # guard
        hit_b = step >= dist_b
        s = np.where(hit_b, dist_b, step)
        # accumulate partial pathlength, move
        L[idx, lay] += s
        x[idx] += s * ux[idx]; y[idx] += s * uy[idx]; z[idx] += s * uzz
        # ---- photons that hit a boundary ----
        hb = idx[hit_b]
        if hb.size:
            lay_hb = layer[hb]; uz_hb = uz[hb]
            going_up = uz_hb < 0
            # snap z to boundary to avoid drift
            z[hb] = np.where(going_up, bounds[lay_hb], bounds_ext[lay_hb + 1])
            # surface check: layer 0 going up to z=0
            at_surface = going_up & (lay_hb == 0)
            # --- surface: Fresnel reflect or exit ---
            su = hb[at_surface]
            if su.size:
                cos_i = -uz[su]                      # uz<0, so >0
                sin_i = np.sqrt(np.clip(1 - cos_i**2, 0, 1))
                sin_t = n_tissue * sin_i
                tir = sin_t >= 1.0
                cos_t = np.sqrt(np.clip(1 - np.minimum(sin_t, 1.0)**2, 0, 1))
                # Fresnel unpolarized reflectance (tissue n_tissue -> air 1)
                rs = (n_tissue*cos_i - 1.0*cos_t) / (n_tissue*cos_i + 1.0*cos_t + 1e-12)
                rp = (n_tissue*cos_t - 1.0*cos_i) / (n_tissue*cos_t + 1.0*cos_i + 1e-12)
                R = 0.5 * (rs**2 + rp**2)
                R = np.where(tir, 1.0, R)
                u = rng.random(su.size)
                reflect = u < R
                # reflected -> flip uz, stay in tissue
                uz[su[reflect]] = -uz[su[reflect]]
                # transmitted -> detected, record exit radius
                ex = su[~reflect]
                detected[ex] = True; active[ex] = False
                exit_r[ex] = np.sqrt(x[ex]**2 + y[ex]**2)
            # --- internal layer boundaries: move to adjacent layer (no scatter) ---
            other = hb[~at_surface]
            if other.size:
                lay_o = layer[other]; up_o = uz[other] < 0
                layer[other] = np.where(up_o, lay_o - 1, lay_o + 1)
                # photons reaching the deep cap (z_max) -> lost
                lost = (~up_o) & (layer[other] >= nlay)
                if lost.any():
                    active[other[lost]] = False
                    layer[other] = np.clip(layer[other], 0, nlay - 1)
        # ---- photons that completed a full step: scatter isotropically (g=0) ----
        sc = idx[~hit_b]
        if sc.size:
            cz = 2*rng.random(sc.size) - 1.0
            phi = 2*np.pi*rng.random(sc.size)
            st = np.sqrt(np.clip(1 - cz**2, 0, 1))
            ux[sc] = st*np.cos(phi); uy[sc] = st*np.sin(phi); uz[sc] = cz
    # ---- score detected photons into SDS annuli ----
    w = np.exp(-(L * mua[None, :]).sum(axis=1))   # Beer-Lambert survival weight
    Lcortex = L[:, cortex_z_layer(bounds, cortex_z):].sum(axis=1)  # layers at/above cortex_z
    Ltot = L.sum(axis=1)
    out = {}
    det = detected & (exit_r >= 0)
    for sds in sds_centers:
        m = det & (np.abs(exit_r - sds) <= half_width)
        ww = w[m]
        out[str(sds)] = dict(
            sum_wLcortex=float((ww*Lcortex[m]).sum()),
            sum_wLtot=float((ww*Ltot[m]).sum()),
            sum_w=float(ww.sum()),
            count=int(m.sum()),
        )
    out['_meta'] = dict(N=N, seed=seed, detected=int(det.sum()),
                        mean_Ltot_all=float(Ltot[det].mean()) if det.any() else 0.0)
    return out

def cortex_z_layer(bounds, cortex_z):
    # index of the first layer whose lower boundary >= cortex_z (i.e., the cortex layer)
    bounds = np.asarray(bounds)
    idxs = np.where(bounds >= cortex_z - 1e-9)[0]
    return int(idxs[0]) if idxs.size else len(bounds)-1

GEOM = {
  '2L': dict(bounds=[0.0, 12.0],        mua=[0.015, 0.020],        musp=[1.0, 0.80],        cortex_z=12.0),
  '3L': dict(bounds=[0.0, 10.0, 12.0],  mua=[0.015, 0.004, 0.020], musp=[1.0, 0.25, 0.80],  cortex_z=12.0),
}

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description=(
            "White Monte Carlo of the layered head slab (Section on the CSF effect). "
            "Scores partial pathlengths to estimate the cortical sensitivity fraction "
            "f_cortex for the two-layer ('2L') or three-layer/CSF ('3L') geometry at "
            "SDS = 25, 30, 35, 40 mm. The ratio gamma = f_cortex(3L)/f_cortex(2L) is the "
            "CSF light-piping factor used in the real-data correction."),
        epilog=(
            "Examples:\n"
            "  python mc_csf.py                       # 3L geometry, paper defaults\n"
            "  python mc_csf.py --geo 2L              # two-layer reference\n"
            "  python mc_csf.py --geo 3L --N 900000 --out csf_3L.json\n"
            "\nf_cortex is stochastic; increase --N to reduce Monte-Carlo noise. The "
            "per-layer optical properties are defined in GEOM (adult-head values)."),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-g', '--geo', choices=['2L', '3L'], default='3L',
                        help="tissue geometry: '2L' (no CSF) or '3L' (with CSF) (default: 3L)")
    parser.add_argument('-N', '--N', type=lambda v: int(float(v)), default=300000,
                        help="number of photons to launch (default: 300000)")
    parser.add_argument('-s', '--seed', type=int, default=1,
                        help="RNG seed (default: 1; converged f_cortex is seed-insensitive)")
    parser.add_argument('-o', '--out', default='mc_csf_fcortex.json',
                        help="output JSON file for the per-SDS results "
                             "(default: mc_csf_fcortex.json)")
    args = parser.parse_args()
    geo, N, seed, outpath = args.geo, args.N, args.seed, args.out
    sds = [25, 30, 35, 40]
    t0=time.time()
    res = run(GEOM[geo], N, seed, sds)
    res['_meta']['geo']=geo; res['_meta']['secs']=round(time.time()-t0,1)
    json.dump(res, open(outpath,'w'), indent=1)
    print(f"wrote results to {outpath}")
    print(f"{geo} N={N} seed={seed} detected={res['_meta']['detected']} "
          f"meanLtot={res['_meta']['mean_Ltot_all']:.1f}mm  ({res['_meta']['secs']}s)")
    for s in sds:
        d=res[str(s)]
        fc = d['sum_wLcortex']/d['sum_wLtot'] if d['sum_wLtot']>0 else float('nan')
        print(f"  SDS={s}mm: count={d['count']:>6d}  f_cortex={fc:.4f}")
