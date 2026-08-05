#!/usr/bin/env python3
"""
Single source of truth for cortical sensitivity fractions.

Runs ONE large anisotropic (g=0.9) Monte-Carlo simulation per (geometry,
wavelength), recording per-photon partial pathlengths, and from that single run
derives BOTH:

  (a) a convergence sweep -- f_cortex re-scored at photon-pathlength cutoffs
      L_max in {400, 500, 800, 1200} mm and detector-annulus half-widths
      {2.0, 2.5, 3.0} mm (plus a separate z_max check) -- to demonstrate the
      long-SDS plateau; and
  (b) batch-based uncertainty -- the detected photons are split into K
      INDEPENDENT batches and f_cortex (a ratio estimator) is recomputed on each
      batch, so the reported mean +/- SD and percentile 95% interval reflect the
      variance of the RATIO across independent batches (not a 3-seed bootstrap).

Wavelength-specific two-layer fractions (760 and 850 nm) and wavelength-specific
CSF light-piping ratios gamma(lambda, SDS) = f_cortex^3L / f_cortex^2L are written
to ONE versioned JSON + CSV.  Both the synthetic-validation and the real-data
pipelines read these values (no hard-coded tables, no 760->850 ratio shortcut).

Run (background, ~1-2 h on 2 cores at the production N):
    python mc_production.py --N 2500000 --batches 16 --out fcortex_production
"""
import argparse, json, time
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from mc_2layer import run

SDS = [25.0, 30.0, 35.0, 38.0, 40.0]
LMAX_SWEEP = [400.0, 500.0, 800.0, 1200.0]
ANNULUS_SWEEP = [2.0, 2.5, 3.0]
PROD_LMAX = 1200.0
PROD_ZMAX = 150.0
PROD_HW = 2.5

# Canonical adult-head optical properties (mm^-1). The 3-layer geometry is the
# 2-layer superficial split into scalp (10 mm) + CSF (2 mm) with identical scalp
# and cortex properties, so gamma = f3L/f2L isolates the CSF light-piping effect.
OPT = {
    760: dict(scalp=(0.0128, 1.08), csf=(0.004, 0.25), cortex=(0.0192, 0.82)),
    850: dict(scalp=(0.0165, 0.86), csf=(0.004, 0.24), cortex=(0.0212, 0.76)),
}

def geom2L(wl):
    o = OPT[wl]
    return dict(bounds=[0.0, 12.0], mua=[o['scalp'][0], o['cortex'][0]],
                musp=[o['scalp'][1], o['cortex'][1]], cortex_z=12.0)

def geom3L(wl, csf_mm=2.0):
    # Superficial depth is held at 12 mm (cortex boundary fixed); the CSF layer
    # occupies the deepest csf_mm of it, so csf_mm=2 -> scalp 0-10 + CSF 10-12,
    # and csf_mm=1 -> scalp 0-11 + CSF 11-12.  Lets gamma be probed vs CSF thickness.
    o = OPT[wl]
    return dict(bounds=[0.0, 12.0 - csf_mm, 12.0],
                mua=[o['scalp'][0], o['csf'][0], o['cortex'][0]],
                musp=[o['scalp'][1], o['csf'][1], o['cortex'][1]], cortex_z=12.0)


import os
NPROC = max(1, min(2, (os.cpu_count() or 2)))


def _worker(task):
    """Top-level worker (picklable) for ProcessPoolExecutor: one MC run."""
    key, geom, N, seed, zmax = task
    raw = run(geom, N, seed, SDS, g=0.9, half_width=PROD_HW, z_max=zmax,
              L_max=PROD_LMAX, return_raw=True)
    return key, raw


def ratio(w, Lcx, Lt, mask):
    sw = (w[mask] * Lt[mask]).sum()
    return float((w[mask] * Lcx[mask]).sum() / sw) if sw > 0 else float('nan')


def score(raw, sds, hw, K):
    """Batch statistics for f_cortex at one SDS annulus."""
    er, Lcx, Lt, w = raw['exit_r'], raw['Lcortex'], raw['Ltot'], raw['w']
    m = np.abs(er - sds) <= hw
    idx = np.where(m)[0]
    if idx.size == 0:
        return None
    f_all = ratio(w, Lcx, Lt, m)
    sw = w[m].sum(); sw2 = (w[m] ** 2).sum()
    neff = float(sw * sw / sw2) if sw2 > 0 else 0.0
    dpf = float((w[m] * Lt[m]).sum() / sw / sds) if sw > 0 else float('nan')
    # K independent batches (disjoint photon subsets)
    rng = np.random.default_rng(0)
    perm = rng.permutation(idx)
    batches = np.array_split(perm, K)
    fb = []
    for b in batches:
        if b.size == 0:
            continue
        bm = np.zeros_like(m); bm[b] = True
        fb.append(ratio(w, Lcx, Lt, bm))
    fb = np.array([x for x in fb if x == x])
    return dict(f_cortex=f_all, batch_mean=float(fb.mean()),
                batch_sd=float(fb.std(ddof=1)) if fb.size > 1 else 0.0,
                ci95=[float(np.percentile(fb, 2.5)), float(np.percentile(fb, 97.5))],
                n_batches=int(fb.size), N_eff_absw=neff, DPF=dpf,
                detected=int(idx.size))


def convergence(raw):
    er, Lcx, Lt, w = raw['exit_r'], raw['Lcortex'], raw['Ltot'], raw['w']
    out = {'Lmax': {}, 'annulus': {}}
    for s in (38.0, 40.0):
        a = np.abs(er - s) <= PROD_HW
        out['Lmax'][str(s)] = {str(L): ratio(w, Lcx, Lt, a & (Lt <= L)) for L in LMAX_SWEEP}
    for s in SDS:
        out['annulus'][str(s)] = {str(hw): ratio(w, Lcx, Lt, np.abs(er - s) <= hw)
                                  for hw in ANNULUS_SWEEP}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('-N', type=lambda v: int(float(v)), default=2500000)
    ap.add_argument('--batches', type=int, default=16)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--zmax_check', type=float, default=220.0)
    ap.add_argument('--zmax_check_N', type=lambda v: int(float(v)), default=800000)
    ap.add_argument('--out', default='fcortex_production')
    args = ap.parse_args()
    t0 = time.time()
    two_layer, csf, conv = {}, {}, {}
    # Build all independent simulations and run them across the available cores
    # (2 at a time), so the four main configs + the z_max check overlap.
    tasks = []
    for wl in (760, 850):
        tasks.append((('2L', wl), geom2L(wl), args.N, args.seed, PROD_ZMAX))
        tasks.append((('3L', wl), geom3L(wl, 2.0), args.N, args.seed, PROD_ZMAX))
        # CSF-thickness robustness: a thinner 1 mm CSF layer (reduced N is fine
        # for this secondary ratio), scored against the same 2-layer fractions.
        tasks.append((('3Lthin', wl), geom3L(wl, 1.0), args.zmax_check_N, args.seed, PROD_ZMAX))
    tasks.append((('zc', 760), geom2L(760), args.zmax_check_N, args.seed, args.zmax_check))
    print(f"[launch] {len(tasks)} independent runs across {NPROC} workers "
          f"(N={args.N}, L_max={PROD_LMAX}, z_max={PROD_ZMAX})", flush=True)
    raws = {}
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        for key, raw in ex.map(_worker, tasks):
            raws[key] = raw
            print(f"[done] {key}  detected={raw['detected']}", flush=True)
    zc = raws.pop(('zc', 760))

    for wl in (760, 850):
        two_layer[wl] = {str(s): score(raws[('2L', wl)], s, PROD_HW, args.batches) for s in SDS}
        conv[wl] = convergence(raws[('2L', wl)])
        csf[wl] = {}
        for s in SDS:
            f2 = two_layer[wl][str(s)]['f_cortex']
            f3s = score(raws[('3L', wl)], s, PROD_HW, args.batches)
            f3 = f3s['f_cortex']
            # gamma uncertainty from batch SDs (relative, added in quadrature)
            rs = np.hypot(two_layer[wl][str(s)]['batch_sd'] / f2 if f2 else 0,
                          f3s['batch_sd'] / f3 if f3 else 0)
            csf[wl][str(s)] = dict(f2L=f2, f3L=f3, gamma=float(f3 / f2) if f2 else None,
                                   gamma_sd=float((f3 / f2) * rs) if f2 else None)
    zc_conv = {str(s): ratio(zc['w'], zc['Lcortex'], zc['Ltot'],
                             np.abs(zc['exit_r'] - s) <= PROD_HW) for s in (38.0, 40.0)}

    # CSF-thickness sweep: gamma with a thinner 1 mm CSF layer (vs the nominal
    # 2 mm), scored against the same converged two-layer fractions.
    csf_thickness = {}
    for wl in (760, 850):
        csf_thickness[wl] = {}
        for s in SDS:
            f2 = two_layer[wl][str(s)]['f_cortex']
            f3t = score(raws[('3Lthin', wl)], s, PROD_HW, args.batches)['f_cortex']
            csf_thickness[wl][str(s)] = dict(f3L_1mm=f3t,
                                             gamma_1mm=float(f3t / f2) if f2 else None)

    result = dict(
        _meta=dict(version="1.0", produced_by="mc_production.py",
                   N_per_config=args.N, n_batches=args.batches, seed=args.seed,
                   g=0.9, L_max=PROD_LMAX, z_max=PROD_ZMAX, half_width=PROD_HW,
                   sds=SDS, optical_properties=OPT,
                   uncertainty="mean +/- SD and 2.5/97.5 percentile across "
                               f"{args.batches} independent photon batches; f_cortex "
                               "is a ratio estimator recomputed per batch",
                   N_eff_note="N_eff_absw = (sum w)^2/sum(w^2) is an absorption-weight "
                              "effective count only; it does not capture pathlength or "
                              "numerator-denominator covariance in the ratio (see batch CI)",
                   convergence_note="L_max/annulus sweeps re-scored from the same run; "
                                    "z_max check is a separate run",
                   csf_note="csf gamma at nominal 2 mm; csf_thickness_1mm holds the "
                            "gamma for a thinner 1 mm CSF layer (superficial depth "
                            "fixed at 12 mm) as a light-piping robustness check",
                   secs=None),
        two_layer=two_layer, csf=csf, csf_thickness_1mm=csf_thickness,
        convergence=conv, zmax_check=zc_conv)
    result['_meta']['secs'] = round(time.time() - t0, 1)
    json.dump(result, open(f"{args.out}.json", 'w'), indent=1)
    with open(f"{args.out}.csv", 'w') as f:
        f.write("geometry,wavelength_nm,SDS_mm,f_cortex,batch_sd,ci95_lo,ci95_hi,"
                "n_batches,N_eff_absw,DPF,detected,gamma,gamma_sd\n")
        for wl in (760, 850):
            for s in SDS:
                d = two_layer[wl][str(s)]; c = csf[wl][str(s)]
                f.write(f"2L,{wl},{s:g},{d['f_cortex']:.5f},{d['batch_sd']:.5f},"
                        f"{d['ci95'][0]:.5f},{d['ci95'][1]:.5f},{d['n_batches']},"
                        f"{d['N_eff_absw']:.1f},{d['DPF']:.3f},{d['detected']},"
                        f"{c['gamma']:.4f},{c['gamma_sd']:.4f}\n")
    # console summary
    print(f"\n=== convergence (2L 760nm, f_cortex vs L_max) ===")
    for s in ('38.0', '40.0'):
        row = conv[760]['Lmax'][s]
        print(f"  SDS={s}: " + "  ".join(f"L{int(float(L))}={row[L]:.4f}" for L in map(str, LMAX_SWEEP)))
    print(f"  z_max {PROD_ZMAX}->{args.zmax_check} (38/40): "
          f"{zc_conv['38.0']:.4f}/{zc_conv['40.0']:.4f}")
    for wl in (760, 850):
        print(f"\n=== two-layer {wl} nm (f_cortex +/- batch SD [95% CI], N_eff, DPF) ===")
        for s in SDS:
            d = two_layer[wl][str(s)]
            print(f"  {s:g}mm: {d['f_cortex']:.4f} +/- {d['batch_sd']:.4f} "
                  f"[{d['ci95'][0]:.4f},{d['ci95'][1]:.4f}]  Neff={d['N_eff_absw']:.0f} DPF={d['DPF']:.2f}")
        print(f"  gamma(CSF,2mm) {wl}: " + "  ".join(
            f"{s:g}:{csf[wl][str(s)]['gamma']:.2f}" for s in SDS))
        print(f"  gamma(CSF,1mm) {wl}: " + "  ".join(
            f"{s:g}:{csf_thickness[wl][str(s)]['gamma_1mm']:.2f}" for s in SDS))
    print(f"\nwrote {args.out}.json / .csv  ({result['_meta']['secs']}s)")


if __name__ == '__main__':
    main()
