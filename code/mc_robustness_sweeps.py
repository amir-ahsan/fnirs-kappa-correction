#!/usr/bin/env python3
"""
Executable SECONDARY robustness sweeps for the fNIRS partial-volume study.

These are the two secondary sensitivity analyses reported in the manuscript
(Figure 3 and the optical-property table).  They are NOT the production
baselines (those come from mc_production.py -> fcortex_production.json); they
deliberately perturb the geometry or the cortical optical properties away from
the production configuration, so they are generated and versioned separately.

This script regenerates them FROM FIRST PRINCIPLES with the two-layer white
Monte Carlo (mc_2layer.run, anisotropic Henyey-Greenstein g=0.9), storing the
RAW per-configuration cortical fractions together with full provenance (seed,
photon count, geometry, optical properties, code commit, command, timestamp).
It writes results/robustness_secondary.json, which the synthetic-validation
plotting reads.  Because it is a fresh, independent Monte-Carlo run with its own
photon count and random stream, its baseline differs modestly from the production
value (about 7% at 30 mm / 760 nm: 0.0597 vs 0.0642).  These sweeps are therefore
used as RELATIVE sensitivity experiments (how f_cortex/kappa_PV move with thickness
and optical properties), NOT as substitute production baselines; the production
baselines come only from mc_production.py -> fcortex_production.json.

  (A) Superficial-thickness sweep: f_cortex and kappa_PV = 1/f_cortex at a fixed
      SDS = 30 mm as the superficial layer thickness varies 8..16 mm (cortex
      boundary moves with it).  Independent MC run per thickness.

  (B) Cortical optical-property sweep at SDS = 30 mm, +/-30 % about the nominal
      cortical mu_a and mu_s':
        - mu_a  is swept by CORRELATED REWEIGHTING of a single photon ensemble
          (w(mu_a) = exp[-(L_sup*mu_a_sup + L_cortex*mu_a_cortex)]), so the mu_a
          curve is low-variance;
        - mu_s' is swept by INDEPENDENT MC runs (scattering cannot be reweighted),
          so it carries the usual few-percent MC noise.

Run (fast; single SDS, moderate N):
    python mc_robustness_sweeps.py -N 500000 --seed 1 --out ../results/robustness_secondary.json
"""
import argparse, json, os, sys, platform, subprocess, hashlib
from datetime import datetime, timezone
from concurrent.futures import ProcessPoolExecutor
import numpy as np
from mc_2layer import run

WL = 760
SDS = 30.0
HW = 2.5
G = 0.9
# Converged pathlength cutoff for this SECONDARY 30 mm sweep. The production
# convergence table shows f_cortex is on its plateau for L_max >= 800 mm; using
# 800 here (vs the production 1200) speeds the thin-superficial configs, whose
# long light-piping paths dominate runtime, with a sub-percent effect at 30 mm.
LMAX = 800.0
NPROC = max(1, min(2, (os.cpu_count() or 2)))
# Nominal 760 nm optical properties (mm^-1), identical to mc_production.OPT[760].
SCALP = (0.0128, 1.08)     # (mu_a, mu_s')
CORTEX = (0.0192, 0.82)
SUP_THICK_NOMINAL = 12.0
THICKNESSES = [8.0, 10.0, 12.0, 14.0, 16.0]   # matches the archived _MC_T validation grid
PCT = [-30, -15, 0, 15, 30]


def _git_commit(cli=None):
    if cli:
        return cli
    try:
        root = os.path.dirname(os.path.abspath(__file__))
        h = subprocess.check_output(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                                    stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.call(["git", "-C", root, "diff", "--quiet"],
                                stderr=subprocess.DEVNULL) != 0
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def _annulus_sums(raw, mua_sup, mua_cortex):
    """Recompute (sum_w*Lcortex, sum_w*Ltot, sum_w, count) at the SDS annulus for
    arbitrary absorption, by reweighting the raw partial pathlengths.  L_sup is
    reconstructed as Ltot - Lcortex (two-layer geometry)."""
    er, Lcx, Lt = raw['exit_r'], raw['Lcortex'], raw['Ltot']
    Lsup = np.clip(Lt - Lcx, 0.0, None)
    w = np.exp(-(Lsup * mua_sup + Lcx * mua_cortex))
    m = np.abs(er - SDS) <= HW
    ww = w[m]
    swl = float((ww * Lcx[m]).sum()); swt = float((ww * Lt[m]).sum())
    return swl, swt, float(ww.sum()), int(m.sum())


def _f_from_raw(raw, mua_sup, mua_cortex):
    swl, swt, _, cnt = _annulus_sums(raw, mua_sup, mua_cortex)
    return (swl / swt if swt > 0 else float('nan')), cnt


def geom2L(sup_thick, mua_cortex=CORTEX[0], musp_cortex=CORTEX[1]):
    return dict(bounds=[0.0, sup_thick], mua=[SCALP[0], mua_cortex],
                musp=[SCALP[1], musp_cortex], cortex_z=sup_thick)


def _run_task(task):
    """Top-level picklable worker: one independent MC run, returned with its label."""
    label, geom, N, seed = task
    raw = run(geom, N, seed, [SDS], g=G, half_width=HW, z_max=150.0, L_max=LMAX,
              return_raw=True)
    return label, raw


def run_sweeps(N, seed):
    # Assemble all INDEPENDENT runs (thickness configs, mu_s' configs, and the
    # single nominal ensemble reused for the mu_a reweighting) and execute them
    # across the available cores.
    tasks = [(("thick", T), geom2L(T), N, seed) for T in THICKNESSES]
    tasks += [(("musp", p), geom2L(SUP_THICK_NOMINAL, musp_cortex=CORTEX[1] * (1 + p / 100.0)),
               N, seed) for p in PCT]
    tasks.append((("base", 0), geom2L(SUP_THICK_NOMINAL), N, seed))
    print(f"  {len(tasks)} independent MC runs across {NPROC} workers "
          f"(N={N}, L_max={LMAX})", flush=True)
    raws = {}
    with ProcessPoolExecutor(max_workers=NPROC) as ex:
        for label, raw in ex.map(_run_task, tasks):
            raws[label] = raw
            print(f"  [done] {label}", flush=True)

    thick_rows = []
    for T in THICKNESSES:
        f, cnt = _f_from_raw(raws[("thick", T)], SCALP[0], CORTEX[0])
        thick_rows.append(dict(sup_thickness_mm=T, f_cortex=float(f),
                               kappa_pv=float(1.0 / f) if f > 0 else None,
                               detected_in_annulus=cnt, N=int(N), seed=int(seed)))
    # mu_a: correlated reweighting of the single nominal ensemble.
    base = raws[("base", 0)]
    mua_rows = []
    for p in PCT:
        mua_c = CORTEX[0] * (1 + p / 100.0)
        f, cnt = _f_from_raw(base, SCALP[0], mua_c)
        mua_rows.append(dict(pct=p, cortex_mua=float(mua_c), f_cortex=float(f),
                             kappa_pv=float(1.0 / f) if f > 0 else None,
                             detected_in_annulus=cnt, estimator="correlated_reweight",
                             N=int(N), seed=int(seed)))
    # mu_s': independent runs.
    musp_rows = []
    for p in PCT:
        f, cnt = _f_from_raw(raws[("musp", p)], SCALP[0], CORTEX[0])
        musp_rows.append(dict(pct=p, cortex_musp=float(CORTEX[1] * (1 + p / 100.0)),
                              f_cortex=float(f), kappa_pv=float(1.0 / f) if f > 0 else None,
                              detected_in_annulus=cnt, estimator="independent_run",
                              N=int(N), seed=int(seed)))
    for r in thick_rows:
        print(f"  thickness {r['sup_thickness_mm']:>4.1f} mm: f_cortex={r['f_cortex']:.4f} "
              f"kappa_PV={r['kappa_pv']:5.2f}", flush=True)
    for r in mua_rows:
        print(f"  mu_a {r['pct']:+d}%: f_cortex={r['f_cortex']:.4f} kappa_PV={r['kappa_pv']:6.3f}", flush=True)
    for r in musp_rows:
        print(f"  mu_s' {r['pct']:+d}%: f_cortex={r['f_cortex']:.4f} kappa_PV={r['kappa_pv']:6.3f}", flush=True)
    return thick_rows, mua_rows, musp_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-N', type=lambda v: int(float(v)), default=500000,
                    help="photons per configuration (canonical frozen value 500000)")
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--out', default='../results/robustness_secondary.json')
    ap.add_argument('--git_commit', default=None,
                    help="Git SHA to record (else auto-detected via git rev-parse)")
    ap.add_argument('--analysis_round', default=None,
                    help="human-readable release label, stored separately from git_commit")
    args = ap.parse_args()

    print("=== SECONDARY robustness sweeps (executable, first-principles MC) ===", flush=True)
    thick, mua_rows, musp_rows = run_sweeps(args.N, args.seed)

    import provenance as _prov
    provenance = _prov.provenance(
        "mc_robustness_sweeps.py",
        analysis_round=args.analysis_round, git_commit=args.git_commit,
        extra=dict(
            engine="mc_2layer.run (Henyey-Greenstein white MC, g=0.9)",
            numpy_version=np.__version__,
            N_per_config=int(args.N), seed=int(args.seed), g=G,
            wavelength_nm=WL, sds_mm=SDS, half_width_mm=HW, L_max=LMAX, z_max=150.0,
            nominal_scalp_mua_musp=list(SCALP), nominal_cortex_mua_musp=list(CORTEX),
            nominal_superficial_thickness_mm=SUP_THICK_NOMINAL))

    out = dict(
        description="SECONDARY robustness sweeps (NOT production baselines): raw "
                    "per-configuration two-layer f_cortex from first-principles white "
                    "Monte Carlo, perturbing superficial thickness (independent runs) "
                    "and cortical optical properties (mu_a by correlated reweighting of "
                    "one ensemble; mu_s' by independent runs). Reported as relative "
                    "sensitivities, not production baselines.",
        provenance=provenance,
        thickness_sweep=dict(sds_mm=SDS, wavelength_nm=WL, rows=thick),
        optical_property_sweep=dict(sds_mm=SDS, wavelength_nm=WL,
                                    mua=mua_rows, musp=musp_rows))
    # payload hash over everything except provenance
    payload = {k: out[k] for k in out if k != "provenance"}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    out["provenance"]["data_sha256"] = hashlib.sha256(blob).hexdigest()

    json.dump(out, open(args.out, 'w'), indent=1)
    print(f"\nwrote {args.out}", flush=True)


if __name__ == '__main__':
    main()
