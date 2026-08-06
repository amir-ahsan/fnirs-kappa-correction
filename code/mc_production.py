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
      LAUNCH-DEFINED batches (batch id = launch_idx mod K) and f_cortex (a ratio
      estimator) is recomputed on each batch, so the reported batch SD, standard
      error (SD/sqrt(K)) and t-based 95% confidence interval of the combined
      estimate reflect the variance of the RATIO across independent batches (not a
      3-seed bootstrap).  Because the batches are keyed by launch index, batch b
      is the same launched-photon cohort across geometries, which makes the CSF
      ratio gamma_b = f3L,b/f2L,b and the z_max difference d_b genuinely PAIRED.

Wavelength-specific two-layer fractions (760 and 850 nm) and wavelength-specific
CSF light-piping ratios gamma(lambda, SDS) = f_cortex^3L / f_cortex^2L are written
to ONE versioned JSON + CSV.  Both the synthetic-validation and the real-data
pipelines read these values (no hard-coded tables, no 760->850 ratio shortcut).

Canonical command used for the manuscript tables (background, ~1-1.5 h on 2 cores):
    python mc_production.py -N 2000000 --batches 16 --out fcortex_production
"""
import argparse, json, time, sys, platform, hashlib, subprocess
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from mc_2layer import run

SCHEMA_VERSION = "2.0"   # bumped when the JSON schema changes (SE/CI/paired fields)


def _git_commit(cli_value=None):
    if cli_value:
        return cli_value
    try:
        here = __file__
        import os
        root = os.path.dirname(os.path.abspath(here))
        h = subprocess.check_output(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                                    stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.call(["git", "-C", root, "diff", "--quiet"],
                                stderr=subprocess.DEVNULL) != 0
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def _payload_sha256(result):
    """SHA-256 over the numeric payload only (excludes _meta), for provenance."""
    payload = {k: result[k] for k in result if k != "_meta"}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()

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


def _t975(dof):
    """Two-sided 95% Student-t critical value; scipy if available, else a table."""
    if dof < 1:
        return float('nan')
    try:
        from scipy.stats import t as _t
        return float(_t.ppf(0.975, dof))
    except Exception:
        tbl = {1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571, 6: 2.447,
               7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228, 11: 2.201, 12: 2.179,
               13: 2.160, 14: 2.145, 15: 2.131, 16: 2.120, 20: 2.086, 30: 2.042}
        if dof in tbl:
            return tbl[dof]
        keys = sorted(tbl)
        return tbl[min(keys, key=lambda k: abs(k - dof))]


def ratio(w, Lcx, Lt, mask):
    sw = (w[mask] * Lt[mask]).sum()
    return float((w[mask] * Lcx[mask]).sum() / sw) if sw > 0 else float('nan')


def score(raw, sds, hw, K):
    """Batch statistics for f_cortex at one SDS annulus.

    Batches are LAUNCH-DEFINED: each detected photon is assigned to batch
    (launch_idx mod K) using its original launch index, not its position in the
    compacted detected-photon array.  Because photons are launched in index order
    with a fixed seed, batch b contains the SAME launched-photon cohort in every
    geometry -- so the per-batch two-/three-layer ratios gamma_b = f3L,b/f2L,b are
    genuinely paired (common random numbers) rather than matched by an arbitrary
    detected-array position.  The batches remain disjoint subsets of the detected
    photons, so they are still valid independent batches for the single-geometry
    SE/CI as well.  (Falls back to a fixed random partition if launch_idx is
    absent, e.g. for a legacy raw dict.)"""
    er, Lcx, Lt, w = raw['exit_r'], raw['Lcortex'], raw['Ltot'], raw['w']
    m = np.abs(er - sds) <= hw
    idx = np.where(m)[0]
    if idx.size == 0:
        return None
    f_all = ratio(w, Lcx, Lt, m)
    sw = w[m].sum(); sw2 = (w[m] ** 2).sum()
    neff = float(sw * sw / sw2) if sw2 > 0 else 0.0
    dpf = float((w[m] * Lt[m]).sum() / sw / sds) if sw > 0 else float('nan')
    # K launch-defined batches: group the annulus photons by (launch_idx mod K),
    # so batch b is the same launched-photon cohort across geometries.
    if 'launch_idx' in raw:
        bt = raw['launch_idx'][idx] % K
        batches = [idx[bt == b] for b in range(K)]
    else:                                     # legacy fallback: fixed random split
        rng = np.random.default_rng(0)
        batches = np.array_split(rng.permutation(idx), K)
    # Keep the per-batch estimates ALIGNED to batch id b (0..K-1): an empty or
    # degenerate batch is stored as NaN rather than dropped, so batch b of the
    # 2-layer run pairs with batch b of the 3-layer run (same launch cohort).
    fb_aligned = []
    for b in batches:
        if b.size == 0:
            fb_aligned.append(float('nan')); continue
        bm = np.zeros_like(m); bm[b] = True
        fb_aligned.append(ratio(w, Lcx, Lt, bm))
    fb_aligned = np.array(fb_aligned, float)
    fb = fb_aligned[np.isfinite(fb_aligned)]   # non-empty batches for summary stats
    nb = int(fb.size)
    bsd = float(fb.std(ddof=1)) if nb > 1 else 0.0
    se = bsd / np.sqrt(nb) if nb > 0 else 0.0          # SE of the combined estimate
    tcrit = _t975(nb - 1)
    # 95% CI for the COMBINED production estimate (t interval on the mean), centred
    # on the full-run ratio f_all.  Distinct from the batch-spread percentiles below.
    ci_t = [float(f_all - tcrit * se), float(f_all + tcrit * se)]
    return dict(f_cortex=f_all, batch_mean=float(fb.mean()),
                batch_sd=bsd, se=float(se), t_crit=float(tcrit),
                ci95=ci_t,
                batch_spread=[float(np.percentile(fb, 2.5)), float(np.percentile(fb, 97.5))],
                batch_estimates=[float(x) for x in fb_aligned],   # length-K, aligned by batch id (NaN=empty); for paired stats
                n_batches=nb, N_eff_absw=neff, DPF=dpf,
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
    ap.add_argument('-N', type=lambda v: int(float(v)), default=2000000,
                    help="photons per configuration (manuscript tables use 2000000)")
    ap.add_argument('--batches', type=int, default=16)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--zmax_check', type=float, default=220.0)
    ap.add_argument('--zmax_check_N', type=lambda v: int(float(v)), default=800000)
    ap.add_argument('--out', default='fcortex_production')
    ap.add_argument('--git_commit', default=None,
                    help="record this commit hash in provenance (else auto-detect)")
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
    # z_max convergence check at the SAME photon count and seed as the production
    # 2L 760 run (common random numbers), so the 150-vs-220 mm difference isolates
    # domain depth rather than sample size.
    tasks.append((('zc', 760), geom2L(760), args.N, args.seed, args.zmax_check))
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
            f2s = two_layer[wl][str(s)]
            f2 = f2s['f_cortex']
            f3s = score(raws[('3L', wl)], s, PROD_HW, args.batches)
            f3 = f3s['f_cortex']
            # gamma uncertainty from the distribution of the per-batch RATIOS
            # gamma_b = f3L,b / f2L,b over LAUNCH-DEFINED batches: batch b of the
            # 2L and 3L runs is the same launched-photon cohort (launch_idx mod B,
            # shared seed = common random numbers), so gamma_b is a genuinely paired
            # estimator.  SE = SD(gamma_b)/sqrt(B_paired).  We also report the
            # independent-propagation SE (quadrature of the 2L and 3L standard
            # errors) as a conservative cross-check that does not rely on pairing.
            b2 = np.array(f2s['batch_estimates']); b3 = np.array(f3s['batch_estimates'])
            pair = np.isfinite(b2) & np.isfinite(b3) & (b2 != 0)
            gb = b3[pair] / b2[pair]
            nbp = int(gb.size)
            gamma_se = float(gb.std(ddof=1) / np.sqrt(nbp)) if nbp > 1 else None
            se2, se3 = f2s['se'], f3s['se']
            gamma_se_indep = (float((f3 / f2) * np.hypot(se2 / f2, se3 / f3))
                              if f2 and f3 else None)
            csf[wl][str(s)] = dict(f2L=f2, f3L=f3, gamma=float(f3 / f2) if f2 else None,
                                   gamma_batch_mean=float(gb.mean()) if nbp else None,
                                   gamma_se=gamma_se, n_paired_batches=nbp,
                                   gamma_se_indep=gamma_se_indep)
    # z_max=220 mm check, scored identically (same N, same seed, same batch
    # partition) as the z_max=150 mm production run, so the batches are truly
    # PAIRED.  Report the mean and SE of the per-batch paired difference
    # d_b = f_150,b - f_220,b (not a quadrature of two separate SEs).
    zc_score = {str(s): score(zc, s, PROD_HW, args.batches) for s in (38.0, 40.0)}
    zc_conv = {}
    for s in ('38.0', '40.0'):
        z2 = zc_score[s]
        f150 = two_layer[760][s]
        b150 = np.array(f150['batch_estimates']); b220 = np.array(z2['batch_estimates'])
        pair = np.isfinite(b150) & np.isfinite(b220)
        db = b150[pair] - b220[pair]
        nb = int(db.size)
        d = float(db.mean())
        d_se = float(db.std(ddof=1) / np.sqrt(nb)) if nb > 1 else 0.0
        zc_conv[s] = dict(f_150=f150['f_cortex'], f_220=z2['f_cortex'],
                          f_220_ci95=z2['ci95'], f_150_ci95=f150['ci95'],
                          delta_paired=d, delta_paired_se=d_se, n_paired_batches=int(nb),
                          within_ci=bool(abs(d) <= 1.96 * d_se) if d_se > 0 else True)

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
        _meta=dict(schema_version=SCHEMA_VERSION, data_version=SCHEMA_VERSION,
                   produced_by="mc_production.py",
                   git_commit=_git_commit(getattr(args, 'git_commit', None)),
                   generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                   command="python mc_production.py " + " ".join(sys.argv[1:]),
                   python_version=platform.python_version(),
                   numpy_version=np.__version__,
                   N_per_config=args.N, n_batches=args.batches, seed=args.seed,
                   g=0.9, L_max=PROD_LMAX, z_max=PROD_ZMAX, half_width=PROD_HW,
                   sds=SDS, optical_properties=OPT,
                   uncertainty=f"f_cortex is a ratio estimator recomputed on each of "
                               f"{args.batches} LAUNCH-DEFINED photon batches (batch id = "
                               "launch_idx mod B; disjoint subsets of the detected photons). "
                               "Reported per SDS: batch_sd (spread of batch estimates), "
                               "se=batch_sd/sqrt(B) (standard error of the combined estimate), "
                               "and ci95 = the t interval f_all +/- t_{.975,B-1}*se (a 95% CI "
                               f"for the combined {args.N:.0e}-photon estimate). batch_spread "
                               "holds the 2.5/97.5 percentiles of the individual batch estimates "
                               "(a prediction range, NOT a CI for the combined run). "
                               "batch_estimates is length B, aligned by batch id (NaN=empty "
                               "batch), so batch b pairs across geometries for the gamma/z_max "
                               "paired statistics.",
                   N_eff_note="N_eff_absw = (sum w)^2/sum(w^2) is an absorption-weight "
                              "effective count only; it does not capture pathlength or "
                              "numerator-denominator covariance in the ratio (see se/ci95)",
                   csf_note="csf gamma at nominal 2 mm. gamma_se is SD(gamma_b)/sqrt(B_paired) "
                            "over the per-batch ratios gamma_b=f3L,b/f2L,b, where batch b is the "
                            "SAME launched-photon cohort in the 2L and 3L runs (launch_idx mod B, "
                            "shared seed) -- a genuinely PAIRED estimator, not a quadrature of two "
                            "separate SEs. gamma_se_indep is the independent-propagation SE "
                            "(quadrature of the 2L and 3L standard errors), reported as a "
                            "conservative cross-check that does not assume pairing. "
                            "csf_thickness_1mm holds gamma for a thinner 1 mm CSF layer.",
                   convergence_note="L_max/annulus sweeps re-scored from the same run; the "
                                    "z_max=220 check uses the SAME N, seed and LAUNCH-DEFINED "
                                    "batches (common random numbers) as the z_max=150 run, so "
                                    "zmax_check reports the mean and SE of the per-batch PAIRED "
                                    "difference d_b=f_150,b-f_220,b.",
                   data_sha256=None, secs=None),
        two_layer=two_layer, csf=csf, csf_thickness_1mm=csf_thickness,
        convergence=conv, zmax_check=zc_conv)
    result['_meta']['secs'] = round(time.time() - t0, 1)
    result['_meta']['data_sha256'] = _payload_sha256(result)
    json.dump(result, open(f"{args.out}.json", 'w'), indent=1)
    with open(f"{args.out}.csv", 'w') as f:
        f.write("geometry,wavelength_nm,SDS_mm,f_cortex,batch_sd,se,ci95_lo,ci95_hi,"
                "n_batches,N_eff_absw,DPF,detected,gamma,gamma_se,gamma_se_indep\n")
        for wl in (760, 850):
            for s in SDS:
                d = two_layer[wl][str(s)]; c = csf[wl][str(s)]
                gsi = c['gamma_se_indep'] if c['gamma_se_indep'] is not None else float('nan')
                f.write(f"2L,{wl},{s:g},{d['f_cortex']:.5f},{d['batch_sd']:.5f},"
                        f"{d['se']:.5f},{d['ci95'][0]:.5f},{d['ci95'][1]:.5f},{d['n_batches']},"
                        f"{d['N_eff_absw']:.1f},{d['DPF']:.3f},{d['detected']},"
                        f"{c['gamma']:.4f},{c['gamma_se']:.4f},{gsi:.4f}\n")
    # console summary
    print(f"\n=== convergence (2L 760nm, f_cortex vs L_max) ===")
    for s in ('38.0', '40.0'):
        row = conv[760]['Lmax'][s]
        print(f"  SDS={s}: " + "  ".join(f"L{int(float(L))}={row[L]:.4f}" for L in map(str, LMAX_SWEEP)))
    print(f"  z_max {PROD_ZMAX}->{args.zmax_check} (matched N={args.N}, seed={args.seed}, paired batches):")
    for s in ('38.0', '40.0'):
        z = zc_conv[s]
        print(f"    SDS={s}: f(150)={z['f_150']:.4f} f(220)={z['f_220']:.4f}  "
              f"paired delta_b={z['delta_paired']:+.5f} +/- {z['delta_paired_se']:.5f} (SE of d_b)  "
              f"{'WITHIN' if z['within_ci'] else 'OUTSIDE'} 95%")
    for wl in (760, 850):
        print(f"\n=== two-layer {wl} nm (f_cortex +/- SE [95% t-CI], batchSD, N_eff, DPF) ===")
        for s in SDS:
            d = two_layer[wl][str(s)]
            print(f"  {s:g}mm: {d['f_cortex']:.4f} +/- {d['se']:.4f} "
                  f"[{d['ci95'][0]:.4f},{d['ci95'][1]:.4f}]  sd={d['batch_sd']:.4f} "
                  f"Neff={d['N_eff_absw']:.0f} DPF={d['DPF']:.2f}")
        print(f"  gamma(CSF,2mm) {wl}: " + "  ".join(
            f"{s:g}:{csf[wl][str(s)]['gamma']:.2f}" for s in SDS))
        print(f"  gamma(CSF,1mm) {wl}: " + "  ".join(
            f"{s:g}:{csf_thickness[wl][str(s)]['gamma_1mm']:.2f}" for s in SDS))
    print(f"\nwrote {args.out}.json / .csv  ({result['_meta']['secs']}s)")


if __name__ == '__main__':
    main()
