"""
In Vivo Proof-of-Concept Demonstration: kappa-Correction of fNIRS Data
======================================================================
Simulates a realistic finger-tapping fNIRS experiment to demonstrate
the kappa-correction pipeline on physiologically realistic data.

The paradigm is INSPIRED BY the BIDS-NIRS-Tapping dataset (Luke et al.,
2021): it parameterizes selected acquisition features (regular 60 blocks of
5 s task / 5 s rest, sampling rate 7.8125 Hz, wavelengths 760 and 850 nm,
SDS ~ 38 mm). It is an idealized block design and does NOT reproduce the real
experiment's irregular inter-onset intervals or its control condition.
This parallel design enables a qualitative comparison between the simulated
pipeline demonstration and the subsequent real-data analysis.

Based on: Sagara N and Ahsan A, "Partial-volume correction in continuous-wave
fNIRS: operating regime, a Monte-Carlo convergence caveat, and a reproducible
pipeline" (submitted to Neurophotonics, 2026)

Dependencies: numpy, matplotlib (no scipy required)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from math import factorial
import os, sys
# Make the single-source loader importable regardless of the working directory.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fcortex_source as _fcs   # reads fcortex_production.json (single source of truth)

# ======================================================================
# Physical Constants and Parameters
# ======================================================================

# Decadic (base-10) molar extinction coefficients in cm^-1 / mM  (Prahl 1999)
# Matched to BIDS-NIRS-Tapping wavelengths: 760 nm and 850 nm
EPS_HBO2 = {760: 0.586, 850: 1.058}
EPS_HBR  = {760: 1.549, 850: 0.691}

# Convert to mm^-1 / uM:  multiply by 1e-4
#   cm^-1/mM = 0.1 mm^-1 / (1000 uM) = 1e-4 mm^-1/uM
UNIT_FACTOR = 1e-4

# DPF values (Scholkmann & Wolf 2013 general equation, age 25)
# Matched to BIDS-NIRS-Tapping wavelengths: 760 nm and 850 nm
DPF = {760: 6.15, 850: 5.089}

# Source-detector separations
# Matched to BIDS-NIRS-Tapping median long-channel SDS (~38 mm)
SDS_LONG  = 38.0   # mm
SDS_SHORT = 10.0   # mm

# Cortical sensitivity fractions: read from the single source of truth
# (fcortex_production.json via fcortex_source.py) -- no hard-coded table here.
# This demo uses the converged two-layer fractions at the long-channel SDS.
F_CORTEX = {wl: _fcs.f_cortex_2L(SDS_LONG, wl) for wl in (760, 850)}
F_CORTEX_SHORT = 0.05  # short channel: almost entirely superficial (illustrative)

# Sampling and timing — inspired by selected BIDS-NIRS-Tapping acquisition parameters
# (a regular block design, NOT a temporal reproduction of the real irregular paradigm)
FS = 7.8125         # Hz (inspired by the BIDS-NIRS-Tapping acquisition rate)
DURATION = 660.0    # seconds (60 blocks x 5s on + 5s off, plus padding)
N_BLOCKS = 60
BLOCK_ON = 5.0      # seconds (inspired by the BIDS-NIRS-Tapping task duration)
BLOCK_OFF = 5.0     # seconds (inspired by the BIDS-NIRS-Tapping rest duration)

# Hemodynamic amplitudes (motor cortex, finger tapping)
HBO2_PEAK = 2.5     # uM
HBR_PEAK  = -0.8    # uM


# ======================================================================
# Helper Functions
# ======================================================================

def hrf(t):
    """Double-gamma hemodynamic response function.
    Parameters: a1=6, b1=1, a2=16, b2=1, c=0.1 (SPM canonical)."""
    a1, b1 = 6, 1
    a2, b2 = 16, 1
    c = 0.1
    h = (t ** (a1 - 1) * np.exp(-t / b1) / (b1 ** a1 * factorial(a1 - 1))
       - c * t ** (a2 - 1) * np.exp(-t / b2) / (b2 ** a2 * factorial(a2 - 1)))
    return h


def bandpass_fft(x, fs, flo=0.01, fhi=0.1):
    """Zero-phase FFT bandpass filter (brick-wall)."""
    X = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(len(x), 1.0 / fs)
    mask = (freqs >= flo) & (freqs <= fhi)
    X[~mask] = 0.0
    return np.fft.irfft(X, n=len(x))


def concentration_to_od(delta_hbo2, delta_hbr, wavelength, sds, dpf):
    """Convert concentration changes (uM) to delta-OD.

    delta_OD(lam) = [eps_HbO2(lam)*dHbO2 + eps_HbR(lam)*dHbR] * L * DPF(lam)
    with eps in mm^-1/uM, dHb in uM, L in mm  -->  dimensionless OD
    """
    eps_hbo2 = EPS_HBO2[wavelength] * UNIT_FACTOR   # mm^-1 / uM
    eps_hbr  = EPS_HBR[wavelength]  * UNIT_FACTOR
    return (eps_hbo2 * delta_hbo2 + eps_hbr * delta_hbr) * sds * dpf


def mbll_inversion(od_760, od_850):
    """Invert MBLL to recover concentrations (uM) from OD.

    b(lam) = OD(lam) / (L * DPF(lam))   in mm^-1
    E * c = b  where E is 2x2 (mm^-1/uM), c is (uM)
    """
    b760 = od_760 / (SDS_LONG * DPF[760])
    b850 = od_850 / (SDS_LONG * DPF[850])

    E = np.array([
        [EPS_HBO2[760] * UNIT_FACTOR, EPS_HBR[760] * UNIT_FACTOR],
        [EPS_HBO2[850] * UNIT_FACTOR, EPS_HBR[850] * UNIT_FACTOR]
    ])
    E_inv = np.linalg.inv(E)

    hbo2 = E_inv[0, 0] * b760 + E_inv[0, 1] * b850
    hbr  = E_inv[1, 0] * b760 + E_inv[1, 1] * b850
    return hbo2, hbr


# ======================================================================
# Simulation
# ======================================================================

def run_simulation():
    np.random.seed(42)
    n_samples = int(DURATION * FS)
    t = np.arange(n_samples) / FS

    # --- Step 1: Create block-design stimulus ---
    # Inspired by selected BIDS-NIRS-Tapping acquisition parameters: 60 regular blocks
    # of 5 s task / 5 s rest (not a temporal reproduction of the real irregular paradigm)
    print("[1/8] Creating block-design stimulus (finger tapping)...")
    stimulus = np.zeros(n_samples)
    pad = 30.0  # initial rest (buffer before first block)
    for b in range(N_BLOCKS):
        t_on  = pad + b * (BLOCK_ON + BLOCK_OFF)
        t_off = t_on + BLOCK_ON
        stimulus[(t >= t_on) & (t < t_off)] = 1.0

    # --- Step 2: Cortical hemodynamic response ---
    print("[2/8] Generating cortical hemodynamic response (HRF convolution)...")
    t_hrf = np.arange(0, 30, 1.0 / FS)
    h = hrf(t_hrf)
    h[0] = 0  # avoid NaN at t=0
    h /= np.max(h)

    hrf_conv = np.convolve(stimulus, h, mode='full')[:n_samples]
    hrf_conv /= np.max(np.abs(hrf_conv)) if np.max(np.abs(hrf_conv)) > 0 else 1.0

    hbo2_true = HBO2_PEAK * hrf_conv   # uM, cortex only
    hbr_true  = HBR_PEAK  * hrf_conv

    # --- Step 3: Superficial physiological noise ---
    # Realistic multi-component systemic noise with slow amplitude modulation.
    # The noise model is tuned to produce R^2_SS values comparable to those
    # observed in the BIDS-NIRS-Tapping dataset (R^2 ~ 0.4-0.7).
    print("[3/8] Simulating superficial physiological noise...")

    # Core systemic oscillations with slow amplitude modulation
    am_cardiac = 1 + 0.3 * np.sin(2*np.pi*0.013*t)
    am_resp    = 1 + 0.3 * np.sin(2*np.pi*0.007*t)
    am_mayer   = 1 + 0.4 * np.sin(2*np.pi*0.005*t)

    cardiac = 0.15 * am_cardiac * np.sin(2*np.pi*1.0*t  + np.random.uniform(0,2*np.pi))
    respir  = 0.25 * am_resp    * np.sin(2*np.pi*0.25*t + np.random.uniform(0,2*np.pi))
    # In-band (0.01-0.1 Hz) systemic amplitudes. This single-channel demo uses
    # a favourable (moderate-contamination) regime so the pipeline stages are
    # clearly visualised; cardiac/respiration are out-of-band and removed by the
    # bandpass. Robust recovery at the true f_cortex requires the channel- and
    # block-averaging used in the real-data analysis (Sec. 4.11).
    mayer   = 0.07 * am_mayer   * np.sin(2*np.pi*0.10*t + np.random.uniform(0,2*np.pi))

    # Additional low-frequency autonomic fluctuations (VLF, 0.02-0.05 Hz)
    vlf = 0.035 * np.sin(2*np.pi*0.04*t + np.random.uniform(0,2*np.pi))
    physio_noise = cardiac + respir + mayer + vlf  # uM equivalent

    hbo2_sup =  physio_noise
    hbr_sup  = -0.3 * physio_noise  # HbR is anti-correlated, smaller

    # --- Step 4: Two-layer forward model (generate measured OD) ---
    print("[4/8] Forward model: two-layer tissue -> measured OD...")
    sigma_od = 2e-4  # instrumental noise std

    # Long channel OD at each wavelength
    od_long = {}
    for lam in [760, 850]:
        od_cortex = concentration_to_od(hbo2_true, hbr_true, lam, SDS_LONG, DPF[lam])
        od_sup    = concentration_to_od(hbo2_sup,  hbr_sup,  lam, SDS_LONG, DPF[lam])
        od_long[lam] = (F_CORTEX[lam] * od_cortex
                        + (1 - F_CORTEX[lam]) * od_sup
                        + sigma_od * np.random.randn(n_samples))

    # Short channel OD (almost purely superficial)
    od_short = {}
    for lam in [760, 850]:
        od_cortex_s = concentration_to_od(hbo2_true, hbr_true, lam, SDS_SHORT, DPF[lam])
        od_sup_s    = concentration_to_od(hbo2_sup,  hbr_sup,  lam, SDS_SHORT, DPF[lam])
        od_short[lam] = (F_CORTEX_SHORT * od_cortex_s
                         + (1 - F_CORTEX_SHORT) * od_sup_s
                         + sigma_od * np.random.randn(n_samples))

    # --- Step 5: Bandpass filter ---
    print("[5/8] Bandpass filtering (0.01 - 0.1 Hz)...")
    for lam in [760, 850]:
        od_long[lam]  = bandpass_fft(od_long[lam],  FS)
        od_short[lam] = bandpass_fft(od_short[lam], FS)

    # --- Step 6: Uncorrected MBLL (baseline) ---
    print("[6/8] MBLL inversion (uncorrected)...")
    hbo2_uncorr, hbr_uncorr = mbll_inversion(od_long[760], od_long[850])

    # --- Step 7: kappa-correction pipeline ---
    print("[7/8] Applying kappa-correction pipeline...")

    # 7a. Short-separation regression (per-wavelength)
    od_ssr = {}
    r2_ss = {}
    for lam in [760, 850]:
        x = od_short[lam]
        y = od_long[lam]
        denom = np.dot(x, x)
        beta = np.dot(x, y) / denom if denom > 1e-20 else 0.0
        od_ssr[lam] = y - beta * x
        # R^2_SS: fraction of variance removed
        var_orig = np.var(y)
        var_resid = np.var(od_ssr[lam])
        r2_ss[lam] = 1 - var_resid / var_orig if var_orig > 1e-20 else 0.0

    # 7b. V_SSR is a DIAGNOSTIC, NOT applied (see KAPPA_PV_ATLAS_CHECK.md
    #     Rec. 2 and the manuscript): SSR already removed the superficial
    #     component; the cortical residual is recovered by kappa_PV alone.
    #     Applying V_SSR here would double-count the dilution.
    od_ssr_corr = {}
    kappa_ssr_vals = {}
    for lam in [760, 850]:
        a_ssr = max(1 - r2_ss[lam], 0.02)     # clamp to avoid blow-up
        kappa_ssr_vals[lam] = 1.0 / a_ssr      # reported as diagnostic only
        od_ssr_corr[lam] = od_ssr[lam]         # V_SSR NOT applied

    # 7c. Partial volume correction: divide by f_cortex (the applied PV factor)
    od_corrected = {}
    for lam in [760, 850]:
        od_corrected[lam] = od_ssr_corr[lam] / F_CORTEX[lam]

    # MBLL on corrected ODs (SSR denoise + PV; V_SSR diagnostic only)
    hbo2_corr, hbr_corr = mbll_inversion(od_corrected[760], od_corrected[850])

    # --- Step 8: Bandpass the ground truth for fair comparison ---
    print("[8/8] Preparing ground truth for comparison...")
    hbo2_gt = bandpass_fft(hbo2_true, FS)
    hbr_gt  = bandpass_fft(hbr_true, FS)

    print("Done!\n")

    # ---- Compute summary statistics ----
    kappa_pv_mean = 0.5 * (1/F_CORTEX[760] + 1/F_CORTEX[850])

    rmse_uncorr_hbo2 = np.sqrt(np.mean((hbo2_uncorr - hbo2_gt)**2))
    rmse_uncorr_hbr  = np.sqrt(np.mean((hbr_uncorr  - hbr_gt)**2))
    rmse_corr_hbo2   = np.sqrt(np.mean((hbo2_corr   - hbo2_gt)**2))
    rmse_corr_hbr    = np.sqrt(np.mean((hbr_corr    - hbr_gt)**2))

    pct_hbo2 = (1 - rmse_corr_hbo2 / rmse_uncorr_hbo2) * 100
    pct_hbr  = (1 - rmse_corr_hbr  / rmse_uncorr_hbr ) * 100

    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"  SDS = {SDS_LONG} mm,  f_cortex(760) = {F_CORTEX[760]},  f_cortex(850) = {F_CORTEX[850]}")
    print(f"  kappa_PV (mean, applied) = {kappa_pv_mean:.3f}")
    print(f"  R^2_SS(760) = {r2_ss[760]:.3f},  R^2_SS(850) = {r2_ss[850]:.3f}")
    print(f"  V_SSR(760) = {kappa_ssr_vals[760]:.3f},  V_SSR(850) = {kappa_ssr_vals[850]:.3f}"
          f"  [per-wavelength diagnostic, NOT applied; no mean/total formed]")
    print()
    print(f"  Ground truth peaks:  HbO2 = {np.max(hbo2_gt):+.3f} uM,  HbR = {np.min(hbr_gt):+.3f} uM")
    print()
    print("  Uncorrected MBLL:")
    print(f"    HbO2 peak = {np.max(hbo2_uncorr):+.3f} uM,  HbR peak = {np.min(hbr_uncorr):+.3f} uM")
    print(f"    RMSE(HbO2) = {rmse_uncorr_hbo2:.4f} uM,  RMSE(HbR) = {rmse_uncorr_hbr:.4f} uM")
    print()
    print("  kappa-Corrected (SSR + PV + MBLL):")
    print(f"    HbO2 peak = {np.max(hbo2_corr):+.3f} uM,  HbR peak = {np.min(hbr_corr):+.3f} uM")
    print(f"    RMSE(HbO2) = {rmse_corr_hbo2:.4f} uM,  RMSE(HbR) = {rmse_corr_hbr:.4f} uM")
    print()
    print(f"  RMSE reduction:  HbO2 = {pct_hbo2:.1f}%,  HbR = {pct_hbr:.1f}%")
    print("=" * 60)

    return dict(
        t=t, stimulus=stimulus,
        hbo2_gt=hbo2_gt, hbr_gt=hbr_gt,
        hbo2_uncorr=hbo2_uncorr, hbr_uncorr=hbr_uncorr,
        hbo2_corr=hbo2_corr, hbr_corr=hbr_corr,
        r2_ss=r2_ss,
        rmse_uncorr_hbo2=rmse_uncorr_hbo2, rmse_uncorr_hbr=rmse_uncorr_hbr,
        rmse_corr_hbo2=rmse_corr_hbo2, rmse_corr_hbr=rmse_corr_hbr,
        pct_hbo2=pct_hbo2, pct_hbr=pct_hbr,
    )


# ======================================================================
# Visualization
# ======================================================================

def make_figures(r, outdir):
    t = r['t']
    stim = r['stimulus']

    # --- Common style ---
    plt.rcParams.update({
        'font.size': 11, 'axes.labelsize': 12,
        'axes.titlesize': 13, 'legend.fontsize': 9,
        'figure.dpi': 150, 'savefig.dpi': 300,
    })

    # ---- Figure 1: Time-series comparison ----
    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    # Shade task blocks
    for ax in axes:
        for i in range(len(stim) - 1):
            if stim[i] > 0.5:
                ax.axvspan(t[i], t[i+1], color='#f0f0f0', lw=0)

    # HbO2
    ax = axes[0]
    ax.plot(t, r['hbo2_gt'],      color='gray',    lw=2.0, label='Ground truth')
    ax.plot(t, r['hbo2_uncorr'],  color='#e74c3c', lw=1.2, ls='--', alpha=0.8, label='Uncorrected MBLL')
    ax.plot(t, r['hbo2_corr'],    color='#c0392b', lw=1.5, label='κ-Corrected')
    ax.set_ylabel('Δ[HbO₂] (μM)')
    ax.set_title('Oxygenated Hemoglobin')
    ax.legend(loc='upper right')
    ax.set_xlim(t[0], t[-1])

    # HbR
    ax = axes[1]
    ax.plot(t, r['hbr_gt'],      color='gray',    lw=2.0, label='Ground truth')
    ax.plot(t, r['hbr_uncorr'],  color='#3498db', lw=1.2, ls='--', alpha=0.8, label='Uncorrected MBLL')
    ax.plot(t, r['hbr_corr'],    color='#2471a3', lw=1.5, label='κ-Corrected')
    ax.set_ylabel('Δ[HbR] (μM)')
    ax.set_xlabel('Time (s)')
    ax.set_title('Deoxygenated Hemoglobin')
    ax.legend(loc='lower right')

    fig.suptitle('Figure 1: Time Series — Uncorrected vs κ-Corrected fNIRS\n'
                 f'(SDS = {SDS_LONG:.0f} mm, finger-tapping block design)',
                 fontsize=14, fontweight='bold', y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path1 = os.path.join(outdir, 'Figure_InVivo_TimeSeries.png')
    fig.savefig(path1)
    plt.close(fig)
    print(f"  Saved: {path1}")

    # ---- Figure 2: Block-averaged HRF ----
    # Matched to real-data analysis: 5 s pre-stimulus, 20 s post-stimulus
    pad_sec = 30.0
    pre_samples = int(5 * FS)  # 5s before onset (same as real-data analysis)
    post_samples = int(20 * FS)  # 20s post-onset (same as real-data analysis)
    epoch_len = pre_samples + post_samples

    def extract_epochs(signal):
        epochs = []
        for b in range(N_BLOCKS):
            onset = int((pad_sec + b * (BLOCK_ON + BLOCK_OFF)) * FS)
            start = onset - pre_samples
            end = start + epoch_len
            if start >= 0 and end <= len(signal):
                ep = signal[start:end].copy()
                ep -= np.mean(ep[:pre_samples])  # baseline correct
                epochs.append(ep)
        return np.array(epochs)

    t_epoch = np.arange(epoch_len) / FS - 5.0  # relative to onset

    ep_gt_hbo2   = extract_epochs(r['hbo2_gt'])
    ep_gt_hbr    = extract_epochs(r['hbr_gt'])
    ep_uc_hbo2   = extract_epochs(r['hbo2_uncorr'])
    ep_uc_hbr    = extract_epochs(r['hbr_uncorr'])
    ep_cr_hbo2   = extract_epochs(r['hbo2_corr'])
    ep_cr_hbr    = extract_epochs(r['hbr_corr'])

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # HbO2 epoch average
    ax = axes[0]
    ax.axvspan(0, BLOCK_ON, color='#f5f5f5', label='Task')
    ax.axhline(0, color='gray', lw=0.5, ls=':')
    ax.plot(t_epoch, ep_gt_hbo2.mean(0),   color='gray',    lw=2.5, label='Ground truth')
    ax.fill_between(t_epoch,
                    ep_gt_hbo2.mean(0) - ep_gt_hbo2.std(0),
                    ep_gt_hbo2.mean(0) + ep_gt_hbo2.std(0),
                    color='gray', alpha=0.15)
    ax.plot(t_epoch, ep_uc_hbo2.mean(0),   color='#e74c3c', lw=1.3, ls='--', label='Uncorrected')
    ax.plot(t_epoch, ep_cr_hbo2.mean(0),   color='#c0392b', lw=1.8, label='κ-Corrected')
    ax.set_xlabel('Time relative to onset (s)')
    ax.set_ylabel('Δ[HbO₂] (μM)')
    ax.set_title('HbO₂ Block-Averaged Response')
    ax.legend(fontsize=8)

    # HbR epoch average
    ax = axes[1]
    ax.axvspan(0, BLOCK_ON, color='#f5f5f5', label='Task')
    ax.axhline(0, color='gray', lw=0.5, ls=':')
    ax.plot(t_epoch, ep_gt_hbr.mean(0),   color='gray',    lw=2.5, label='Ground truth')
    ax.fill_between(t_epoch,
                    ep_gt_hbr.mean(0) - ep_gt_hbr.std(0),
                    ep_gt_hbr.mean(0) + ep_gt_hbr.std(0),
                    color='gray', alpha=0.15)
    ax.plot(t_epoch, ep_uc_hbr.mean(0),   color='#3498db', lw=1.3, ls='--', label='Uncorrected')
    ax.plot(t_epoch, ep_cr_hbr.mean(0),   color='#2471a3', lw=1.8, label='κ-Corrected')
    ax.set_xlabel('Time relative to onset (s)')
    ax.set_ylabel('Δ[HbR] (μM)')
    ax.set_title('HbR Block-Averaged Response')
    ax.legend(fontsize=8)

    fig.suptitle('Figure 2: Block-Averaged Hemodynamic Response\n'
                 f'(mean ± SD across {N_BLOCKS} blocks)',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    path2 = os.path.join(outdir, 'Figure_InVivo_HRF.png')
    fig.savefig(path2, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path2}")

    # ---- Figure 3: Summary bar chart ----
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    # Peak amplitudes
    ax = axes[0]
    gt_peaks  = [np.max(r['hbo2_gt']),  abs(np.min(r['hbr_gt']))]
    uc_peaks  = [np.max(r['hbo2_uncorr']), abs(np.min(r['hbr_uncorr']))]
    cr_peaks  = [np.max(r['hbo2_corr']),   abs(np.min(r['hbr_corr']))]
    x = np.arange(2)
    w = 0.25
    ax.bar(x - w, gt_peaks,  w, color='gray',    label='Ground truth', edgecolor='black', lw=0.5)
    ax.bar(x,     uc_peaks,  w, color='#e8a0a0',  label='Uncorrected',  edgecolor='black', lw=0.5)
    ax.bar(x + w, cr_peaks,  w, color='#2ecc71',  label='κ-Corrected',  edgecolor='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(['|Δ[HbO₂]| peak', '|Δ[HbR]| peak'])
    ax.set_ylabel('Peak amplitude (μM)')
    ax.set_title('Peak Concentration Recovery')
    ax.legend()

    # RMSE
    ax = axes[1]
    rmse_vals = [
        [r['rmse_uncorr_hbo2'], r['rmse_uncorr_hbr']],
        [r['rmse_corr_hbo2'],   r['rmse_corr_hbr']],
    ]
    ax.bar(x - 0.15, rmse_vals[0], 0.3, color='#e8a0a0',  label='Uncorrected', edgecolor='black', lw=0.5)
    ax.bar(x + 0.15, rmse_vals[1], 0.3, color='#2ecc71',  label='κ-Corrected',  edgecolor='black', lw=0.5)
    ax.set_xticks(x)
    ax.set_xticklabels(['HbO₂', 'HbR'])
    ax.set_ylabel('RMSE (μM)')
    ax.set_title('Root-Mean-Square Error')
    ax.legend()

    # Add improvement % annotation
    for i, (uc, cr) in enumerate(zip(rmse_vals[0], rmse_vals[1])):
        pct = (1 - cr / uc) * 100
        ax.annotate(f'{pct:.1f}% ↓',
                    xy=(i + 0.15, cr), xytext=(i + 0.4, cr + 0.05),
                    fontsize=9, fontweight='bold', color='#27ae60',
                    arrowprops=dict(arrowstyle='->', color='#27ae60'))

    fig.suptitle('Figure 3: Quantitative Improvement with κ-Correction\n'
                 f'(SDS = {SDS_LONG:.0f} mm, simulated finger-tapping task)',
                 fontsize=14, fontweight='bold', y=1.02)
    fig.tight_layout()
    path3 = os.path.join(outdir, 'Figure_InVivo_Summary.png')
    fig.savefig(path3, bbox_inches='tight')
    plt.close(fig)
    print(f"  Saved: {path3}")


# ======================================================================
# Main
# ======================================================================

if __name__ == '__main__':
    print("=" * 60)
    print("In Vivo Proof-of-Concept Demonstration")
    print("kappa-Correction of fNIRS Data")
    print("=" * 60)
    print()

    results = run_simulation()

    outdir = os.path.dirname(os.path.abspath(__file__))
    print("\nGenerating figures...")
    make_figures(results, outdir)
    print("\nAll done!")
