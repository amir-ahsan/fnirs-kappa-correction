#!/usr/bin/env python3
"""
LEGACY / NOT USED FOR THE MANUSCRIPT TABLES.

This earlier uncertainty/provenance driver has been SUPERSEDED by
mc_production.py, which is the single source of truth for every f_cortex and CSF
ratio reported in the manuscript (written to fcortex_production.json with schema
version, provenance, launch-defined paired batches, and a payload hash).
mc_uncertainty.py is retained only for historical reference; it is not imported
by any analysis script and must not be used to regenerate the manuscript numbers.

Monte-Carlo uncertainty / provenance driver for the cortical sensitivity
fraction f_cortex and the CSF light-piping ratio gamma.

It uses the SAME anisotropic Henyey-Greenstein transport (g=0.9) as mc_2layer.py
for BOTH wavelengths and for BOTH the two-layer and three-layer (CSF) geometries,
so the two-layer table and the CSF ratio are computed under identical transport
assumptions.

For every (geometry, wavelength, SDS) it reports, aggregated over several
independent seeds:
  * f_cortex        = sum(w*L_cortex) / sum(w*L_total)         (per-seed, then mean +/- SD)
  * 95% CI          = across-seed bootstrap percentile interval
  * detected count  = photons in the SDS annulus (mean over seeds)
  * N_eff           = (sum w)^2 / sum(w^2)   effective (weighted) sample size
  * DPF             = (sum(w*L_total)/sum(w)) / SDS   weighted mean pathlength / SDS
The CSF ratio gamma = f_cortex(3L)/f_cortex(2L_matched) is computed at several
CSF thicknesses with propagated uncertainty.  Raw per-seed accumulators and the
aggregated table are written to JSON and CSV.

Note on "photons per separation": each seed launches one ensemble and bins exit
radii into SDS annuli; the per-SDS detected and N_eff counts below make the
effective sample size at each separation explicit (it is far smaller than the
launched count, especially at long SDS).
"""
import argparse, json, time
import numpy as np
from mc_2layer import run  # anisotropic HG (g=0.9) layered-slab engine

# ---- optical properties (adult head; paper Table of optical properties) ----
# Two-layer "sensitivity-table" geometries (reproduce the adopted f_cortex table)
GEOM_2L = {
    760: dict(bounds=[0.0, 12.0], mua=[0.0128, 0.0192], musp=[1.08, 0.82], cortex_z=12.0),
    850: dict(bounds=[0.0, 12.0], mua=[0.0165, 0.0212], musp=[0.86, 0.76], cortex_z=12.0),
}
# CSF-ratio geometries: matched 2L vs 3L differing ONLY by an explicit CSF layer,
# both with g=0.9, so gamma isolates the CSF light-piping effect.
def geom_matched(wl, csf_mm):
    scalp = dict(mua=0.015, musp=1.0) if wl == 760 else dict(mua=0.019, musp=0.90)
    csf   = dict(mua=0.004, musp=0.25)
    ctx   = dict(mua=0.020, musp=0.80) if wl == 760 else dict(mua=0.024, musp=0.72)
    scalp_thick = 12.0 - csf_mm
    g2 = dict(bounds=[0.0, 12.0], mua=[scalp['mua'], ctx['mua']],
              musp=[scalp['musp'], ctx['musp']], cortex_z=12.0)
    g3 = dict(bounds=[0.0, scalp_thick, 12.0],
              mua=[scalp['mua'], csf['mua'], ctx['mua']],
              musp=[scalp['musp'], csf['musp'], ctx['musp']], cortex_z=12.0)
    return g2, g3

SDS = [25, 30, 35, 38, 40]


def aggregate(per_seed, sds_list):
    """per_seed: list of run() outputs. Returns per-SDS aggregated stats."""
    out = {}
    for s in sds_list:
        k = str(s)
        fcs, dpfs, neffs, counts = [], [], [], []
        for r in per_seed:
            d = r[k]
            if d['sum_wLtot'] > 0 and d['sum_w'] > 0:
                fcs.append(d['sum_wLcortex'] / d['sum_wLtot'])
                dpfs.append((d['sum_wLtot'] / d['sum_w']) / s)
                neffs.append(d['sum_w']**2 / d['sum_w2'] if d.get('sum_w2', 0) > 0 else 0.0)
                counts.append(d['count'])
        fcs = np.array(fcs)
        if fcs.size == 0:
            out[k] = None; continue
        # across-seed bootstrap percentile CI (resample the seeds)
        rng = np.random.default_rng(12345)
        boot = [np.mean(rng.choice(fcs, fcs.size, replace=True)) for _ in range(2000)] \
            if fcs.size > 1 else [fcs[0]]
        out[k] = dict(
            f_cortex=float(fcs.mean()), f_cortex_sd=float(fcs.std(ddof=1)) if fcs.size > 1 else 0.0,
            ci95=[float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))],
            kappa_PV=float(1.0 / fcs.mean()),
            DPF=float(np.mean(dpfs)), N_eff=float(np.mean(neffs)),
            detected=float(np.mean(counts)), n_seeds=int(fcs.size))
    return out


def run_multi(geom, N, seeds, g=0.9):
    return [run(geom, N, sd, SDS, g=g) for sd in seeds]


def main():
    ap = argparse.ArgumentParser(description="Multi-seed Monte-Carlo uncertainty driver for f_cortex and CSF gamma.")
    ap.add_argument('-N', type=lambda v: int(float(v)), default=200000, help="photons per seed (default 2e5)")
    ap.add_argument('--seeds', type=int, default=4, help="number of independent seeds (default 4)")
    ap.add_argument('--csf', type=float, nargs='+', default=[1.0, 2.0, 3.0], help="CSF thicknesses (mm)")
    ap.add_argument('--out', default='mc_uncertainty', help="output basename (.json/.csv)")
    args = ap.parse_args()
    seeds = list(range(1, args.seeds + 1))
    t0 = time.time()
    result = {'_meta': dict(N=args.N, seeds=seeds, g=0.9, sds=SDS)}

    # 1) two-layer sensitivity table, both wavelengths
    twolayer = {}
    for wl in (760, 850):
        agg = aggregate(run_multi(GEOM_2L[wl], args.N, seeds), SDS)
        twolayer[wl] = agg
        print(f"[2L {wl}nm] " + "  ".join(
            f"{s}:{agg[str(s)]['f_cortex']:.4f}±{agg[str(s)]['f_cortex_sd']:.4f}"
            f"(Neff{agg[str(s)]['N_eff']:.0f},DPF{agg[str(s)]['DPF']:.1f})" for s in SDS))
    result['two_layer'] = twolayer

    # 2) CSF ratio gamma at several thicknesses (760 nm; 850 at 2 mm)
    csf = {}
    for wl in (760, 850):
        for t in (args.csf if wl == 760 else [2.0]):
            g2, g3 = geom_matched(wl, t)
            a2 = aggregate(run_multi(g2, args.N, seeds), SDS)
            a3 = aggregate(run_multi(g3, args.N, seeds), SDS)
            gam = {}
            for s in SDS:
                k = str(s)
                f2, f3 = a2[k]['f_cortex'], a3[k]['f_cortex']
                # propagate: relative sd added in quadrature
                rs = np.hypot(a2[k]['f_cortex_sd'] / f2 if f2 else 0,
                              a3[k]['f_cortex_sd'] / f3 if f3 else 0)
                gam[k] = dict(gamma=float(f3 / f2) if f2 else None,
                              gamma_sd=float((f3 / f2) * rs) if f2 else None,
                              f2=f2, f3=f3)
            csf[f"{wl}nm_csf{t:g}mm"] = gam
            print(f"[CSF {wl}nm t={t:g}mm] " + "  ".join(
                f"{s}:g={gam[str(s)]['gamma']:.2f}±{gam[str(s)]['gamma_sd']:.2f}" for s in SDS))
    result['csf_gamma'] = csf
    result['_meta']['secs'] = round(time.time() - t0, 1)

    json.dump(result, open(f"{args.out}.json", 'w'), indent=1)
    # flat CSV of the two-layer table
    with open(f"{args.out}.csv", 'w') as f:
        f.write("geometry,wavelength_nm,SDS_mm,f_cortex,f_cortex_sd,ci95_lo,ci95_hi,kappa_PV,DPF,N_eff,detected,n_seeds\n")
        for wl in (760, 850):
            for s in SDS:
                d = twolayer[wl][str(s)]
                f.write(f"2L,{wl},{s},{d['f_cortex']:.5f},{d['f_cortex_sd']:.5f},"
                        f"{d['ci95'][0]:.5f},{d['ci95'][1]:.5f},{d['kappa_PV']:.3f},"
                        f"{d['DPF']:.3f},{d['N_eff']:.1f},{d['detected']:.0f},{d['n_seeds']}\n")
    print(f"\nwrote {args.out}.json and {args.out}.csv  ({result['_meta']['secs']}s)")


if __name__ == '__main__':
    main()
