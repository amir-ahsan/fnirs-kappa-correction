#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Group-level kappa-correction analysis on all five BIDS-NIRS-Tapping subjects.

Reuses the per-subject pipeline functions defined in fnirs_kappa_realdata_analysis.py
(here imported as `realdata_helpers`), loops over sub-01 through sub-05, and
aggregates kappa factors and corrected/uncorrected peak amplitudes into
group-level statistics.

Outputs:
    group_summary.json      -- machine-readable summary
    group_summary.csv       -- per-subject table
    group_summary.txt       -- human-readable summary
    group_block_average.png -- HbO2 and HbR block-averaged HRFs (mean +- SE
                                across subjects)
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import fnirs_kappa_realdata_analysis as rh

SCRIPT_DIR = Path(__file__).resolve().parent
DATA_DIR = SCRIPT_DIR / "BIDS-NIRS-Tapping-data"
OUT_DIR = SCRIPT_DIR / "outputs"
OUT_DIR.mkdir(exist_ok=True)

SUBJECTS = ["01", "02", "03", "04", "05"]


# ---------------------------------------------------------------------------
# Per-subject pipeline (extracted from fnirs_kappa_realdata_analysis.main)
# ---------------------------------------------------------------------------
def analyze_subject(subject: str):
    """Run the kappa-correction pipeline on one subject and return results."""
    print(f"\n{'='*70}\n  Subject sub-{subject}\n{'='*70}")
    raw = rh.load_subject(DATA_DIR, subject)
    fs = raw.info["sfreq"]

    long_picks, short_picks, wavelengths = rh.split_channels(raw)
    if len(wavelengths) < 2:
        raise RuntimeError("Need at least two wavelengths")
    wl1, wl2 = wavelengths[0], wavelengths[1]

    # Fill in missing constants for any unusual wavelengths
    for wl in (wl1, wl2):
        if wl not in rh.EXTINCTION_COEFFICIENTS:
            nearest = min(rh.EXTINCTION_COEFFICIENTS, key=lambda x: abs(x - wl))
            rh.EXTINCTION_COEFFICIENTS[wl] = rh.EXTINCTION_COEFFICIENTS[nearest]
        if wl not in rh.DPF:
            nearest = min(rh.DPF, key=lambda x: abs(x - wl))
            rh.DPF[wl] = rh.DPF[nearest]

    # SDS
    dists = rh.source_detector_distances(raw.info)
    long_all = np.concatenate([long_picks[wl1], long_picks[wl2]])
    sds_mm = float(np.median(dists[long_all])) * 1000.0

    # f_cortex
    f_cortex = rh.get_f_cortex(sds_mm)
    for wl in (wl1, wl2):
        if wl not in f_cortex:
            nearest = min(f_cortex, key=lambda x: abs(x - wl))
            f_cortex[wl] = f_cortex[nearest]

    # Pipeline
    raw_od = rh.optical_density(raw)
    od_data = raw_od.get_data().T
    od_filt = rh.bandpass_filter(od_data, fs, lo=0.01, hi=0.1)

    # SSR + V_SSR
    od_ssr = od_filt.copy()
    r2_by_wl, kappa_ssr_by_wl = {}, {}
    for wl in (wl1, wl2):
        lp = long_picks[wl]
        sp = short_picks.get(wl, np.array([], dtype=int))
        if len(sp) == 0:
            r2_by_wl[wl] = np.zeros(len(lp))
            kappa_ssr_by_wl[wl] = 1.0
            continue
        corrected, betas, r2 = rh.ssr_per_wavelength(od_filt[:, lp], od_filt[:, sp])
        od_ssr[:, lp] = corrected
        r2_by_wl[wl] = r2
        kappa_ssr_by_wl[wl] = float(np.mean(rh.compute_kappa_ssr(r2)))

    # V_SSR is a DIAGNOSTIC, NOT applied.  1/(1-R^2_SS) is an inverse
    # residual-VARIANCE ratio, not an amplitude-restoration factor: the
    # residual-to-original standard-deviation (RMS-fluctuation) ratio would be
    # sqrt(1-R^2), but neither quantity identifies cortical task-signal loss.
    # kappa_PV = 1/f_cortex corrects only the optical dilution of the cortical
    # component that survives SSR; it does not reconstruct signal SSR removed.
    od_ssr_corrected = od_ssr.copy()

    # PV correction
    od_pv = od_ssr_corrected.copy()
    kappa_pvs = {}
    for wl in (wl1, wl2):
        f = f_cortex[wl]
        kpv = 1.0 / f
        kappa_pvs[wl] = kpv
        lp = long_picks[wl]
        od_pv[:, lp] = od_ssr_corrected[:, lp] / f

    kappa_pv_mean = float(np.mean(list(kappa_pvs.values())))
    # V_SSR is a per-wavelength diagnostic (not applied); no cross-wavelength
    # mean or kappa_total is formed (the 760/850 nm values differ too much for
    # their average to be meaningful).

    # MBLL inversion (channel-by-channel)
    def _label(name):
        return name.split()[0] if " " in name else name
    sd_wl1 = {_label(raw.ch_names[i]): i for i in long_picks[wl1]}
    sd_wl2 = {_label(raw.ch_names[i]): i for i in long_picks[wl2]}
    common = sorted(set(sd_wl1) & set(sd_wl2))

    # od_pv_only: partial-volume correction applied WITHOUT the SSR step, so the
    # net scaling factor can be decomposed into an SSR-only contribution and a
    # PV-only contribution (referee point 3.2).  od_ssr (= od_ssr_corrected) is
    # the SSR-only path (no PV); od_pv is the full SSR+PV path.
    od_pv_only = od_filt.copy()
    for wl in (wl1, wl2):
        od_pv_only[:, long_picks[wl]] = od_filt[:, long_picks[wl]] / f_cortex[wl]

    hbo2_u_all, hbr_u_all, hbo2_c_all, hbr_c_all = [], [], [], []
    hbo2_ssr_all, hbo2_pv_all = [], []
    for label in common:
        i1, i2 = sd_wl1[label], sd_wl2[label]
        hbo2_u, hbr_u = rh.mbll_inversion(od_filt[:, i1], od_filt[:, i2], wl1, wl2, sds_mm)
        hbo2_c, hbr_c = rh.mbll_inversion(od_pv[:, i1], od_pv[:, i2], wl1, wl2, sds_mm)
        hbo2_s, _ = rh.mbll_inversion(od_ssr_corrected[:, i1], od_ssr_corrected[:, i2], wl1, wl2, sds_mm)
        hbo2_p, _ = rh.mbll_inversion(od_pv_only[:, i1], od_pv_only[:, i2], wl1, wl2, sds_mm)
        hbo2_u_all.append(hbo2_u); hbr_u_all.append(hbr_u)
        hbo2_c_all.append(hbo2_c); hbr_c_all.append(hbr_c)
        hbo2_ssr_all.append(hbo2_s); hbo2_pv_all.append(hbo2_p)

    hbo2_u = np.mean(hbo2_u_all, axis=0); hbr_u = np.mean(hbr_u_all, axis=0)
    hbo2_c = np.mean(hbo2_c_all, axis=0); hbr_c = np.mean(hbr_c_all, axis=0)
    hbo2_ssr = np.mean(hbo2_ssr_all, axis=0); hbo2_pv = np.mean(hbo2_pv_all, axis=0)

    # Block average around tapping events
    events = rh.get_tapping_events(raw)
    samples = events[:, 0]
    ba = {}
    for name, sig in [("hbo2_u", hbo2_u), ("hbr_u", hbr_u), ("hbo2_c", hbo2_c),
                      ("hbr_c", hbr_c), ("hbo2_ssr", hbo2_ssr), ("hbo2_pv", hbo2_pv)]:
        t_ba, m, s, n = rh.block_average(sig, samples, fs, pre=5.0, post=20.0)
        ba["time"] = t_ba
        ba[name + "_mean"] = m
        ba[name + "_std"] = s
        ba["n_blocks"] = n

    # Peaks are taken within the physiological response window
    # (rh.PEAK_WINDOW_S) rather than over the whole epoch, matching the
    # single-subject script.
    hbo2_peak_u = rh.windowed_peak_abs(ba["hbo2_u_mean"], ba["time"])
    hbo2_peak_c = rh.windowed_peak_abs(ba["hbo2_c_mean"], ba["time"])
    hbr_peak_u = rh.windowed_peak_abs(ba["hbr_u_mean"], ba["time"])
    hbr_peak_c = rh.windowed_peak_abs(ba["hbr_c_mean"], ba["time"])
    # SSR-only and PV-only HbO2 peaks: decompose the net scaling into the SSR
    # denoising factor (raw -> SSR-only) and the PV amplification factor
    # (SSR-only -> full).  ssr_scale * pv_scale ~= net hbo2_scale.
    hbo2_peak_ssr = rh.windowed_peak_abs(ba["hbo2_ssr_mean"], ba["time"])
    hbo2_peak_pv = rh.windowed_peak_abs(ba["hbo2_pv_mean"], ba["time"])

    return {
        "subject": subject,
        "sds_mm": sds_mm,
        "wavelengths": [wl1, wl2],
        "n_channels": len(common),
        "n_blocks": int(ba["n_blocks"]),
        "f_cortex": {str(k): float(v) for k, v in f_cortex.items() if k in (wl1, wl2)},
        "kappa_pv": {str(k): float(v) for k, v in kappa_pvs.items()},
        "kappa_pv_mean": kappa_pv_mean,
        "kappa_ssr": {str(k): float(v) for k, v in kappa_ssr_by_wl.items()},
        "kappa_ssr_760": float(kappa_ssr_by_wl.get(760, float("nan"))),
        "kappa_ssr_850": float(kappa_ssr_by_wl.get(850, float("nan"))),
        "hbo2_peak_uncorr_uM": hbo2_peak_u,
        "hbo2_peak_ssr_only_uM": hbo2_peak_ssr,
        "hbo2_peak_pv_only_uM": hbo2_peak_pv,
        "hbo2_peak_corr_uM": hbo2_peak_c,
        "hbo2_ssr_scale": hbo2_peak_ssr / hbo2_peak_u if hbo2_peak_u > 0 else float("nan"),
        "hbo2_pv_scale": hbo2_peak_c / hbo2_peak_ssr if hbo2_peak_ssr > 0 else float("nan"),
        "hbo2_scale": hbo2_peak_c / hbo2_peak_u if hbo2_peak_u > 0 else float("nan"),
        "hbr_peak_uncorr_uM": hbr_peak_u,
        "hbr_peak_corr_uM": hbr_peak_c,
        "hbr_scale": hbr_peak_c / hbr_peak_u if hbr_peak_u > 0 else float("nan"),
        "_ba": ba,  # keep for figure
    }


# ---------------------------------------------------------------------------
# Group aggregation
# ---------------------------------------------------------------------------
def aggregate(results):
    """Compute mean +/- SD across subjects for the headline metrics."""
    keys = ["kappa_pv_mean", "kappa_ssr_760", "kappa_ssr_850",
            "hbo2_peak_uncorr_uM", "hbo2_peak_ssr_only_uM", "hbo2_peak_pv_only_uM",
            "hbo2_peak_corr_uM", "hbo2_ssr_scale", "hbo2_pv_scale", "hbo2_scale",
            "hbr_peak_uncorr_uM", "hbr_peak_corr_uM", "hbr_scale"]
    summary = {}
    for k in keys:
        vals = np.array([r[k] for r in results], dtype=float)
        summary[k] = {"mean": float(np.mean(vals)), "sd": float(np.std(vals, ddof=1)),
                      "n": int(len(vals)), "values": vals.tolist()}
    return summary


def write_outputs(per_subject, group):
    # CSV
    cols = ["subject", "sds_mm", "n_channels", "n_blocks", "kappa_pv_mean",
            "kappa_ssr_760", "kappa_ssr_850",
            "hbo2_peak_uncorr_uM", "hbo2_peak_ssr_only_uM", "hbo2_peak_pv_only_uM",
            "hbo2_peak_corr_uM", "hbo2_ssr_scale", "hbo2_pv_scale", "hbo2_scale",
            "hbr_peak_uncorr_uM", "hbr_peak_corr_uM", "hbr_scale"]
    csv_lines = [",".join(cols)]
    for r in per_subject:
        csv_lines.append(",".join(f"{r[c]:.4f}" if isinstance(r[c], float) else str(r[c]) for c in cols))
    (OUT_DIR / "group_summary.csv").write_text("\n".join(csv_lines))

    # JSON (drop _ba)
    clean = [{k: v for k, v in r.items() if k != "_ba"} for r in per_subject]
    (OUT_DIR / "group_summary.json").write_text(
        json.dumps({"per_subject": clean, "group": group}, indent=2)
    )

    # Human-readable
    lines = ["Group-level kappa-correction analysis (BIDS-NIRS-Tapping, N=5)",
             "=" * 70, ""]
    lines.append(f"{'Subject':<10}{'SDS(mm)':>10}{'kPV':>8}{'V_SSR760':>9}{'V_SSR850':>9}"
                 f"{'HbO2 unc':>10}{'HbO2 cor':>10}{'scale':>8}"
                 f"{'HbR unc':>10}{'HbR cor':>10}{'scale':>8}")
    lines.append("-" * 110)
    for r in per_subject:
        lines.append(
            f"{'sub-' + r['subject']:<10}{r['sds_mm']:>10.1f}"
            f"{r['kappa_pv_mean']:>8.3f}{r['kappa_ssr_760']:>9.3f}{r['kappa_ssr_850']:>9.3f}"
            f"{r['hbo2_peak_uncorr_uM']:>10.3f}{r['hbo2_peak_corr_uM']:>10.3f}{r['hbo2_scale']:>8.2f}"
            f"{r['hbr_peak_uncorr_uM']:>10.3f}{r['hbr_peak_corr_uM']:>10.3f}{r['hbr_scale']:>8.2f}"
        )
    lines.append("-" * 110)
    lines.append(f"{'Mean +- SD':<10}{'':>10}"
                 f"{group['kappa_pv_mean']['mean']:>8.3f}"
                 f"{group['kappa_ssr_760']['mean']:>9.3f}"
                 f"{group['kappa_ssr_850']['mean']:>9.3f}"
                 f"{group['hbo2_peak_uncorr_uM']['mean']:>10.3f}"
                 f"{group['hbo2_peak_corr_uM']['mean']:>10.3f}"
                 f"{group['hbo2_scale']['mean']:>8.2f}"
                 f"{group['hbr_peak_uncorr_uM']['mean']:>10.3f}"
                 f"{group['hbr_peak_corr_uM']['mean']:>10.3f}"
                 f"{group['hbr_scale']['mean']:>8.2f}")
    lines.append(f"{'(SD)':<10}{'':>10}"
                 f"{group['kappa_pv_mean']['sd']:>8.3f}"
                 f"{group['kappa_ssr_760']['sd']:>9.3f}"
                 f"{group['kappa_ssr_850']['sd']:>9.3f}"
                 f"{group['hbo2_peak_uncorr_uM']['sd']:>10.3f}"
                 f"{group['hbo2_peak_corr_uM']['sd']:>10.3f}"
                 f"{group['hbo2_scale']['sd']:>8.2f}"
                 f"{group['hbr_peak_uncorr_uM']['sd']:>10.3f}"
                 f"{group['hbr_peak_corr_uM']['sd']:>10.3f}"
                 f"{group['hbr_scale']['sd']:>8.2f}")
    (OUT_DIR / "group_summary.txt").write_text("\n".join(lines))
    print("\n".join(lines))


def make_group_figure(per_subject):
    """Plot block-averaged HRFs across subjects."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for r in per_subject:
        ba = r["_ba"]; t = ba["time"]
        axes[0].plot(t, ba["hbo2_u_mean"], color="lightcoral", alpha=0.55, lw=1)
        axes[0].plot(t, ba["hbo2_c_mean"], color="firebrick", alpha=0.85, lw=1.4)
        axes[1].plot(t, ba["hbr_u_mean"], color="lightskyblue", alpha=0.55, lw=1)
        axes[1].plot(t, ba["hbr_c_mean"], color="steelblue", alpha=0.85, lw=1.4)
    # group means
    t = per_subject[0]["_ba"]["time"]
    hbo2_u_grp = np.mean([r["_ba"]["hbo2_u_mean"] for r in per_subject], axis=0)
    hbo2_c_grp = np.mean([r["_ba"]["hbo2_c_mean"] for r in per_subject], axis=0)
    hbr_u_grp = np.mean([r["_ba"]["hbr_u_mean"] for r in per_subject], axis=0)
    hbr_c_grp = np.mean([r["_ba"]["hbr_c_mean"] for r in per_subject], axis=0)
    axes[0].plot(t, hbo2_u_grp, "k--", lw=2, label="Uncorrected (group mean)")
    axes[0].plot(t, hbo2_c_grp, "k-", lw=2.2, label=r"$\kappa$-corrected (group mean)")
    axes[1].plot(t, hbr_u_grp, "k--", lw=2, label="Uncorrected (group mean)")
    axes[1].plot(t, hbr_c_grp, "k-", lw=2.2, label=r"$\kappa$-corrected (group mean)")
    for ax, ttl in zip(axes, [r"$\Delta$[HbO$_2$] ($\mu$M)", r"$\Delta$[HbR] ($\mu$M)"]):
        ax.axvspan(0, 5, color="grey", alpha=0.15, label="Task")
        ax.set_xlabel("Time relative to onset (s)")
        ax.set_ylabel(ttl)
        ax.legend(loc="best", fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("BIDS-NIRS-Tapping group analysis (N = 5; per-subject thin lines, group mean bold)")
    fig.tight_layout()
    fig.savefig(OUT_DIR / "group_block_average.png", dpi=150)
    print(f"  Figure saved -> {OUT_DIR / 'group_block_average.png'}")


def main():
    per_subject = []
    for sub in SUBJECTS:
        try:
            per_subject.append(analyze_subject(sub))
        except Exception as e:
            print(f"  ERROR on sub-{sub}: {e}")
    if not per_subject:
        sys.exit(1)
    group = aggregate(per_subject)
    write_outputs(per_subject, group)
    make_group_figure(per_subject)
    print(f"\nAll outputs saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
