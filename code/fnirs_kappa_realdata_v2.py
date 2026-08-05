#!/usr/bin/env python3
"""
Upgraded in-vivo analysis of BIDS-NIRS-Tapping (addresses reviewer Section 4.4).

Improvements over fnirs_kappa_realdata_analysis.py:
  * per-CHANNEL source-detector separation and per-channel f_cortex (not one median);
  * condition-resolved (Tapping/Left vs Tapping/Right) block averaging with
    CONTRALATERAL channel selection (left hand -> right hemisphere, x>0;
    right hand -> left hemisphere, x<0; convention validated by lateralization);
  * NEAREST short-separation channel as the SSR regressor for each long channel
    (not the global short-channel mean);
  * quality control: motion correction (TDDR) + scalp-coupling-index (SCI) rejection;
  * a fixed pre-registered response-window MEAN amplitude estimator (not a max-|.|,
    which is positively biased after large amplification).

Only kappa_PV = 1/f_cortex(per-channel SDS, CSF-augmented) is applied; kappa_SSR is
reported as a variance-removal diagnostic (not applied).

Requires: mne, mne-bids, numpy, scipy.  Reuses constants/helpers from
fnirs_kappa_realdata_analysis.py.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import mne
from mne_bids import BIDSPath, read_raw_bids
from mne.preprocessing.nirs import (optical_density, source_detector_distances,
                                    scalp_coupling_index,
                                    temporal_derivative_distribution_repair)

import fnirs_kappa_realdata_analysis as rd  # constants + get_f_cortex

mne.set_log_level("ERROR")

# Pre-registered response windows (s after onset) for the window-mean estimator.
HBO_WINDOW = (4.0, 10.0)   # canonical HbO peak band for a 5 s block
HBR_WINDOW = (7.0, 13.0)   # HbR reaches its (negative) extremum later
SCI_THRESHOLD = 0.5        # scalp-coupling-index QC threshold
PRE, POST = 5.0, 17.0      # epoch window (s)


def _base(name):  # "S1_D1 760" -> "S1_D1"
    return name.split(" ")[0]


def load_subject(root: Path, sub: str) -> mne.io.Raw:
    bp = BIDSPath(subject=sub, task="tapping", datatype="nirs", root=root,
                  suffix="nirs", extension=".snirf")
    return read_raw_bids(bp, verbose="ERROR").load_data()


def analyze_subject(root: Path, sub: str) -> dict:
    raw = load_subject(root, sub)
    sf = raw.info["sfreq"]

    # --- optical density, motion correction (TDDR), QC by SCI (pre-filter) ---
    od = optical_density(raw)
    od = temporal_derivative_distribution_repair(od)
    sci = scalp_coupling_index(od)
    good = {name: (sci[i] >= SCI_THRESHOLD) for i, name in enumerate(od.ch_names)}
    od.filter(0.01, 0.1, l_trans_bandwidth=0.005, h_trans_bandwidth=0.02,
              verbose="ERROR")

    dists_m = source_detector_distances(od.info)
    data = od.get_data()  # (n_ch, T), natural-log OD

    # channel bookkeeping
    chans = []
    for i, ch in enumerate(od.info["chs"]):
        name = ch["ch_name"]
        try:
            wl = int(round(ch["loc"][9])) if not np.isnan(ch["loc"][9]) else None
        except Exception:
            wl = None
        if wl not in (760, 850):
            wl = 760 if name.endswith("760") else (850 if name.endswith("850") else None)
        chans.append(dict(idx=i, name=name, base=_base(name), wl=wl,
                          sds_mm=float(dists_m[i]) * 1000.0,
                          x=float(ch["loc"][0]), mid=ch["loc"][:3].copy(),
                          good=bool(good[name])))
    # short channels (SDS < 15 mm) per wavelength, good only
    shorts = {wl: [c for c in chans if c["wl"] == wl and c["sds_mm"] < 15 and c["good"]]
              for wl in (760, 850)}
    # long S-D bases present & good at BOTH wavelengths
    by_base = {}
    for c in chans:
        by_base.setdefault(c["base"], {})[c["wl"]] = c
    long_bases = [b for b, d in by_base.items()
                  if 760 in d and 850 in d and d[760]["sds_mm"] >= 25
                  and d[760]["good"] and d[850]["good"]]

    def nearest_short(c):
        cand = shorts[c["wl"]]
        if not cand:
            return None
        dd = [np.linalg.norm(c["mid"] - s["mid"]) for s in cand]
        return cand[int(np.argmin(dd))]

    # --- per-channel SSR (nearest short), PV, MBLL ---
    E = np.array([[rd.EXTINCTION_COEFFICIENTS[wl]["HbO2"],
                   rd.EXTINCTION_COEFFICIENTS[wl]["HbR"]] for wl in (760, 850)]) * rd.UNIT_CONVERSION
    E_inv = np.linalg.pinv(E)

    chan_out = []  # per long S-D: HbO/HbR time series (uncorr, corr), x, sds
    kssr = {760: [], 850: []}
    for b in long_bases:
        od_ss = {}
        od_un = {}
        fcx = {}
        ok = True
        for wl in (760, 850):
            c = by_base[b][wl]
            y = data[c["idx"]] / rd.LN10           # decadic OD
            s = nearest_short(c)
            if s is None:
                ok = False; break
            x = data[s["idx"]] / rd.LN10
            # OLS SSR (zero-intercept on band-passed, ~zero-mean signals)
            beta = float(np.dot(x, y) / np.dot(x, x)) if np.dot(x, x) > 0 else 0.0
            resid = y - beta * x
            r2 = 1.0 - np.var(resid) / np.var(y) if np.var(y) > 0 else 0.0
            r2 = float(np.clip(r2, 0.0, 0.99))
            kssr[wl].append(1.0 / (1.0 - r2))
            od_un[wl] = y
            od_ss[wl] = resid
            fcx[wl] = rd.get_f_cortex(c["sds_mm"])[wl]   # per-channel, CSF-augmented
        if not ok:
            continue
        sds = by_base[b][760]["sds_mm"]
        Lu = np.array([rd.DPF[wl] * sds for wl in (760, 850)])
        # uncorrected concentrations
        ODu = np.vstack([od_un[760] / Lu[0], od_un[850] / Lu[1]])
        Cu = E_inv @ ODu                      # (2, T): HbO, HbR
        # corrected: SSR + per-channel PV then MBLL
        ODc = np.vstack([(od_ss[760] / fcx[760]) / Lu[0],
                         (od_ss[850] / fcx[850]) / Lu[1]])
        Cc = E_inv @ ODc
        chan_out.append(dict(base=b, x=by_base[b][760]["x"], sds=sds,
                             kpv=1.0 / np.mean([fcx[760], fcx[850]]),
                             f760=fcx[760], f850=fcx[850],
                             kpv760=1.0 / fcx[760], kpv850=1.0 / fcx[850],
                             hbo_u=Cu[0], hbr_u=Cu[1], hbo_c=Cc[0], hbr_c=Cc[1]))

    # --- events ---
    ev, eid = mne.events_from_annotations(raw, verbose="ERROR")
    onsets = {"left": ev[ev[:, 2] == eid["Tapping/Left"], 0],
              "right": ev[ev[:, 2] == eid["Tapping/Right"], 0]}
    pre_s, post_s = int(PRE * sf), int(POST * sf)
    t = np.arange(-pre_s, post_s) / sf

    def block_mean(sig, ons):
        segs = []
        for o in ons:
            lo, hi = o - pre_s, o + post_s
            if lo < 0 or hi > sig.shape[0]:
                continue
            e = sig[lo:hi] - sig[lo:lo + pre_s].mean()
            segs.append(e)
        return np.mean(segs, axis=0) if segs else np.zeros(pre_s + post_s)

    def win_amp(mean_trace, window):
        w0 = pre_s + int(window[0] * sf); w1 = pre_s + int(window[1] * sf)
        return float(mean_trace[w0:w1].mean())

    # contralateral selection: left hand -> right hemi (x>0); right hand -> left hemi (x<0)
    contra = {"left": lambda c: c["x"] > 0, "right": lambda c: c["x"] < 0}
    res = {"hbo_u": [], "hbo_c": [], "hbr_u": [], "hbr_c": [], "kpv": [], "sds": []}
    tr = {"hbo_u": [], "hbo_c": [], "hbr_u": [], "hbr_c": []}   # block-average traces
    for cond in ("left", "right"):
        sel = [c for c in chan_out if contra[cond](c)]
        if not sel:
            continue
        for c in sel:
            m_hbo_u = block_mean(c["hbo_u"], onsets[cond])
            m_hbo_c = block_mean(c["hbo_c"], onsets[cond])
            m_hbr_u = block_mean(c["hbr_u"], onsets[cond])
            m_hbr_c = block_mean(c["hbr_c"], onsets[cond])
            res["hbo_u"].append(win_amp(m_hbo_u, HBO_WINDOW))
            res["hbo_c"].append(win_amp(m_hbo_c, HBO_WINDOW))
            res["hbr_u"].append(win_amp(m_hbr_u, HBR_WINDOW))
            res["hbr_c"].append(win_amp(m_hbr_c, HBR_WINDOW))
            res["kpv"].append(c["kpv"]); res["sds"].append(c["sds"])
            tr["hbo_u"].append(m_hbo_u); tr["hbo_c"].append(m_hbo_c)
            tr["hbr_u"].append(m_hbr_u); tr["hbr_c"].append(m_hbr_c)
    traces = {k: (np.mean(v, axis=0) if v else np.zeros(pre_s + post_s)) for k, v in tr.items()}
    traces["t"] = t

    # Concentrations from the inversion are already in micromolar (E carries the
    # mM^-1 cm^-1 -> mm^-1 uM^-1 conversion), so no further scaling is applied.
    hbo_u = float(np.mean(res["hbo_u"])); hbo_c = float(np.mean(res["hbo_c"]))
    hbr_u = float(np.mean(res["hbr_u"])); hbr_c = float(np.mean(res["hbr_c"]))
    return dict(subject=sub, n_long=len(chan_out),
                n_contra=len(res["hbo_u"]),
                sds_mean=float(np.mean(res["sds"])),
                kpv_mean=float(np.mean(res["kpv"])),
                kssr760=float(np.mean(kssr[760])), kssr850=float(np.mean(kssr[850])),
                hbo_uncorr_uM=hbo_u, hbo_corr_uM=hbo_c,
                hbr_uncorr_uM=hbr_u, hbr_corr_uM=hbr_c,
                hbo_scale=hbo_c / hbo_u if hbo_u else float("nan"),
                channels=[dict(base=c["base"], sds_mm=round(c["sds"], 2),
                               hemisphere=("right" if c["x"] > 0 else "left"),
                               f_cortex_760=round(c["f760"], 4), f_cortex_850=round(c["f850"], 4),
                               kappa_pv_760=round(c["kpv760"], 3), kappa_pv_850=round(c["kpv850"], 3))
                          for c in chan_out],
                traces=traces)


def make_figures(rows, script_dir):
    """Regenerate the real-data figures from the upgraded contralateral pipeline."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    t = rows[0]["traces"]["t"]
    def stack(key):
        return np.vstack([r["traces"][key] for r in rows])
    fig_dir = script_dir  # written next to the script; copied into manuscript/ later
    # --- Figure 5: single-subject (sub-01) contralateral block-averaged HRF ---
    r0 = rows[0]["traces"]
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    ax[0].axvspan(0, 5, color="0.9"); ax[1].axvspan(0, 5, color="0.9")
    ax[0].plot(t, r0["hbo_u"], "--", color="tab:red", label="uncorrected")
    ax[0].plot(t, r0["hbo_c"], "-", color="tab:red", label="κ-corrected")
    ax[0].set_title("Δ[HbO₂] (Subject 01, contralateral)"); ax[0].legend(fontsize=8)
    ax[1].plot(t, r0["hbr_u"], "--", color="tab:blue", label="uncorrected")
    ax[1].plot(t, r0["hbr_c"], "-", color="tab:blue", label="κ-corrected")
    ax[1].set_title("Δ[HbR] (Subject 01, contralateral)"); ax[1].legend(fontsize=8)
    for a in ax:
        a.set_xlabel("time (s)"); a.set_ylabel("Δc (µM)"); a.axhline(0, color="k", lw=0.5)
    fig.tight_layout(); fig.savefig(fig_dir / "figure5_realdata_hrf.png", dpi=150); plt.close(fig)
    # --- Figure 7: group block-average (per-subject light, group-mean bold) ---
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, chrom, col in ((ax[0], "hbo", "tab:red"), (ax[1], "hbr", "tab:blue")):
        a.axvspan(0, 5, color="0.9")
        U, C = stack(chrom + "_u"), stack(chrom + "_c")
        for i in range(U.shape[0]):
            a.plot(t, U[i], color=col, alpha=0.25, lw=0.8)
            a.plot(t, C[i], color=col, alpha=0.25, lw=0.8, ls="--")
        a.plot(t, U.mean(0), color="k", lw=2, label="group uncorrected")
        a.plot(t, C.mean(0), color=col, lw=2.5, label="group κ-corrected")
        a.set_title(f"Δ[{'HbO₂' if chrom=='hbo' else 'HbR'}] group (N=5, contralateral)")
        a.set_xlabel("time (s)"); a.set_ylabel("Δc (µM)")
        a.axhline(0, color="k", lw=0.5); a.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(fig_dir / "figure7_group_block_average.png", dpi=150); plt.close(fig)
    # --- Figure 4: per-subject amplification (uncorrected vs corrected HbO) ---
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    u = np.array([r["hbo_uncorr_uM"] for r in rows]); c = np.array([r["hbo_corr_uM"] for r in rows])
    ax.scatter(u, c, s=60, color="tab:red", zorder=3)
    for r, xu, yc in zip(rows, u, c):
        ax.annotate(f"sub-{r['subject']}", (xu, yc), fontsize=7, xytext=(4, 3), textcoords="offset points")
    kpv = np.mean([r["kpv_mean"] for r in rows])
    xl = np.array([0, max(u) * 1.15])
    ax.plot(xl, kpv * xl, "--", color="0.5", label=f"slope = mean κ_PV = {kpv:.1f}")
    ax.plot(xl, xl, ":", color="0.7", label="identity")
    ax.set_xlabel("uncorrected contralateral Δ[HbO₂] (µM)")
    ax.set_ylabel("κ-corrected Δ[HbO₂] (µM)")
    ax.set_title("Per-subject amplification (N=5, SDS≈38 mm)")
    ax.legend(fontsize=8); ax.set_xlim(0, xl[1]); ax.set_ylim(0, max(c) * 1.15)
    fig.tight_layout(); fig.savefig(fig_dir / "figure4_realdata_timeseries.png", dpi=150); plt.close(fig)
    # --- Figure 6: group summary bars ---
    fig, ax = plt.subplots(figsize=(5.6, 4.4))
    gU = [abs(np.mean([r["hbo_uncorr_uM"] for r in rows])), abs(np.mean([r["hbr_uncorr_uM"] for r in rows]))]
    gC = [abs(np.mean([r["hbo_corr_uM"] for r in rows])), abs(np.mean([r["hbr_corr_uM"] for r in rows]))]
    eU = [np.std([r["hbo_uncorr_uM"] for r in rows], ddof=1), np.std([r["hbr_uncorr_uM"] for r in rows], ddof=1)]
    eC = [np.std([r["hbo_corr_uM"] for r in rows], ddof=1), np.std([r["hbr_corr_uM"] for r in rows], ddof=1)]
    xpos = np.arange(2); w = 0.35
    ax.bar(xpos - w / 2, gU, w, yerr=eU, capsize=4, label="uncorrected", color="0.7")
    ax.bar(xpos + w / 2, gC, w, yerr=eC, capsize=4, label="κ-corrected", color="tab:red")
    ax.set_xticks(xpos); ax.set_xticklabels(["|Δ[HbO₂]|", "|Δ[HbR]|"])
    ax.set_ylabel("group-mean amplitude (µM)")
    ax.set_title("Group summary (N=5, contralateral, window-mean)")
    ax.legend(fontsize=9)
    fig.tight_layout(); fig.savefig(fig_dir / "figure6_realdata_summary.png", dpi=150); plt.close(fig)
    print("wrote figures 4,5,6,7 (real-data, upgraded pipeline)")


def main():
    script_dir = Path(__file__).resolve().parent
    root = script_dir / "BIDS-NIRS-Tapping-data"
    if not (root.exists() and any(root.glob("sub-*"))):
        root = rd.download_dataset(root)
    rows = []
    for sub in ["01", "02", "03", "04", "05"]:
        r = analyze_subject(root, sub)
        rows.append(r)
        print(f"sub-{sub}: n_contra={r['n_contra']:2d} SDS={r['sds_mean']:.1f}mm "
              f"kPV={r['kpv_mean']:.2f} kSSR(760/850)={r['kssr760']:.2f}/{r['kssr850']:.2f}  "
              f"HbO {r['hbo_uncorr_uM']:+.3f}->{r['hbo_corr_uM']:+.3f}uM ({r['hbo_scale']:.2f}x)  "
              f"HbR {r['hbr_uncorr_uM']:+.3f}->{r['hbr_corr_uM']:+.3f}uM")
    arr = lambda k: np.array([r[k] for r in rows])
    print("\nGROUP (contralateral, window-mean, N=5):")
    for k, lab in [("hbo_uncorr_uM", "HbO uncorr"), ("hbo_corr_uM", "HbO corr"),
                   ("hbr_uncorr_uM", "HbR uncorr"), ("hbr_corr_uM", "HbR corr"),
                   ("kpv_mean", "kPV"), ("kssr760", "kSSR760"), ("kssr850", "kSSR850")]:
        print(f"  {lab:11s}: {arr(k).mean():+.3f} +/- {arr(k).std(ddof=1):.3f}")
    out = dict(per_subject=rows,
               group=dict(hbo_uncorr=float(arr("hbo_uncorr_uM").mean()),
                          hbo_uncorr_sd=float(arr("hbo_uncorr_uM").std(ddof=1)),
                          hbo_corr=float(arr("hbo_corr_uM").mean()),
                          hbo_corr_sd=float(arr("hbo_corr_uM").std(ddof=1)),
                          hbr_uncorr=float(arr("hbr_uncorr_uM").mean()),
                          hbr_corr=float(arr("hbr_corr_uM").mean()),
                          kpv=float(arr("kpv_mean").mean()),
                          kssr760=float(arr("kssr760").mean()), kssr760_sd=float(arr("kssr760").std(ddof=1)),
                          kssr850=float(arr("kssr850").mean()), kssr850_sd=float(arr("kssr850").std(ddof=1))))
    make_figures(rows, script_dir)
    for r in out["per_subject"]:
        r.pop("traces", None)
    json.dump(out, open(script_dir / "realdata_v2_summary.json", "w"), indent=1)
    print("wrote realdata_v2_summary.json")


if __name__ == "__main__":
    main()
