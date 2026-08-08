#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fnirs_kappa_realdata_analysis.py
=================================

Real-data demonstration of the kappa-correction framework for functional
near-infrared spectroscopy (fNIRS). This script accompanies the manuscript:

    Sagara N and Ahsan A (2026). "Partial-volume correction in continuous-wave
    fNIRS: operating regime, a Monte-Carlo convergence caveat, and a reproducible
    pipeline." Neurophotonics (submitted).

The pipeline applies four processing stages to optical-density (OD)
data acquired from a publicly available finger-tapping experiment:

    1. Short-separation regression (SSR) at the OD level, performed
       independently for each wavelength.
    2. SSR variance-removal diagnostic: V_SSR(lam) = 1 / (1 - R^2_SS(lam)),
       estimated per wavelength from the regression statistics and REPORTED
       as a diagnostic of the long-channel variance removed by (coupled to)
       the short-channel regression, NOT of superficial contamination
       remaining after SSR. It is NOT
       applied to the data, and no cross-wavelength mean or kappa_total is
       formed: applying V_SSR on top of the (correctly large) kappa_PV
       would double-count the dilution, and the per-wavelength values differ
       too much for their mean to be meaningful.
    3. Partial-volume (PV) correction: OD_cortex = OD_ssr / f_cortex(lam),
       where f_cortex is the cortical sensitivity fraction from the converged
       two-layer Monte Carlo (mc_2layer.py).
    4. Modified Beer-Lambert Law (MBLL) inversion to recover changes in
       oxygenated (HbO2) and deoxygenated (HbR) haemoglobin concentration.

Dataset
-------
BIDS-NIRS-Tapping (Luke et al., 2021)
    Repository : https://github.com/rob-luke/BIDS-NIRS-Tapping
    DOI        : 10.5281/zenodo.5529797
    Paradigm   : Alternating left/right finger-tapping blocks
    Wavelengths: 760 nm, 850 nm
    Subjects   : 5 (sub-01 through sub-05)

The dataset is downloaded automatically on first execution.

Requirements
------------
    Python >= 3.9
    numpy, scipy, matplotlib, mne, mne-nirs, mne-bids

Install with::

    pip install numpy scipy matplotlib mne mne-nirs mne-bids

Usage
-----
Run from the command line::

    python fnirs_kappa_realdata_analysis.py

Three publication-quality figures are saved to the script directory:

    - Figure_RealData_TimeSeries.png
    - Figure_RealData_HRF.png
    - Figure_RealData_Summary.png

License
-------
Released under the MIT License; see the root LICENSE file for the full terms.
(The manuscript text and figures are separately licensed under CC BY 4.0; see
LICENSE-manuscript.md.)

Authors
-------
Neth Sagara, Amir Ahsan
Irvine Valley College
"""

# ── Standard library ────────────────────────────────────────────────────
from __future__ import annotations

import hashlib
import os
import sys
import urllib.request
import warnings
import zipfile
from pathlib import Path
from typing import Any

# ── Third-party ─────────────────────────────────────────────────────────
import numpy as np
from numpy.typing import NDArray
from scipy import signal as sp_signal
import matplotlib.pyplot as plt

import mne
from mne.preprocessing.nirs import optical_density, source_detector_distances
from mne_bids import BIDSPath, read_raw_bids

# ── Suppress verbose MNE output ────────────────────────────────────────
mne.set_log_level("WARNING")
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ========================================================================
#  Constants
# ========================================================================

#: Molar extinction coefficients (cm-1 / mM, decadic / base-10).
#: Source: Prahl / Oregon Medical Laser Center (OMLC) compilation.
EXTINCTION_COEFFICIENTS: dict[int, dict[str, float]] = {
    760: {"HbO2": 0.586, "HbR": 1.549},
    850: {"HbO2": 1.058, "HbR": 0.691},
}

#: Conversion factor from cm-1 / mM to mm-1 / uM.
UNIT_CONVERSION: float = 1e-4

#: Natural-log -> base-10 conversion.  MNE's optical_density() returns
#: natural-log optical density (ln I/I0), whereas the Prahl extinction
#: coefficients above are decadic (base-10).  Dividing the OD by ln(10)
#: converts it to decadic absorbance so the two are consistent; omitting
#: this overestimates every recovered concentration by ln(10) ~ 2.303x.
LN10: float = float(np.log(10.0))

#: Differential pathlength factors (dimensionless).
#: Source: Scholkmann and Wolf (2013), general equation evaluated at age 25.
#: DPF(760) = 6.15, DPF(850) = 5.09 (the earlier 5.35 was inconsistent with the
#: cited equation and with the synthetic pipeline, which uses the formula value).
DPF: dict[int, float] = {760: 6.15, 850: 5.089}

#: Cortical sensitivity fractions come from the SINGLE production Monte-Carlo
#: source (fcortex_production.json via fcortex_source): direct, wavelength-specific
#: two-layer fractions at BOTH 760 and 850 nm and wavelength-specific CSF ratios
#: gamma(lambda, SDS).  There is NO hard-coded 760-nm table and NO fixed 760->850
#: ratio here -- the earlier 1.02 shortcut and single-gamma multiplier are removed,
#: so the in-vivo correction propagates the full revised Monte-Carlo calculation.
import fcortex_source as _fcs

APPLY_CSF_LIGHTPIPING: bool = True


def csf_gamma(sds_mm: float, wavelength_nm: int = 760) -> float:
    """Wavelength-specific CSF light-piping ratio gamma(SDS, lambda) (production MC)."""
    return _fcs.gamma_csf(sds_mm, wavelength_nm)

#: Maximum source-detector distance (m) for a channel to be classified
#: as "short-separation".
SS_THRESHOLD: float = 0.015  # 15 mm

#: Matplotlib styling for publication-quality figures.
_PLOT_PARAMS: dict[str, Any] = {
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "legend.fontsize": 9,
}


# ========================================================================
#  Helper functions
# ========================================================================

def get_f_cortex(sds_mm: float) -> dict[int, float]:
    """Per-wavelength cortical sensitivity fractions for an arbitrary SDS.

    Reads direct 760- and 850-nm two-layer fractions from the single production
    Monte-Carlo source and applies the corresponding wavelength-specific CSF
    light-piping ratio gamma(lambda, SDS) for real heads.  No 760->850 ratio
    shortcut and no single gamma multiplier are used.

    Returns the FULL floating-point interpolated fraction (no rounding): callers
    that apply the quantitative kappa_PV correction must use full precision, and
    any rounding is done only at display/export time. Rounding here would quantise
    every corrected amplitude to the 4th-decimal grid of f_cortex.
    """
    out = {}
    for wl in (760, 850):
        f = (_fcs.f_cortex_invivo(sds_mm, wl) if APPLY_CSF_LIGHTPIPING
             else _fcs.f_cortex_2L(sds_mm, wl))
        out[wl] = float(f)
    return out


# ========================================================================
#  1. Dataset acquisition
# ========================================================================

# PINNED, immutable dataset source. We do NOT download the moving GitHub `master`
# branch (mutable, would break reproducibility). Instead we resolve the frozen,
# citable Zenodo release (a Zenodo record is immutable once published) via the
# Zenodo REST API, download the archive it lists, and VERIFY the Zenodo-published
# checksum (md5). This makes the input data reproducible and self-verifying without
# hard-coding a possibly-stale URL or hash. If the network is unavailable, a clear
# manual-download instruction is printed. DATASET_SHA256, if set, additionally
# enforces a SHA-256 on the downloaded archive.
DATASET_DOI = "10.5281/zenodo.5529797"          # rob-luke/BIDS-NIRS-Tapping (Luke et al., 2021)
DATASET_ZENODO_RECORD = "5529797"
DATASET_SHA256 = None                            # optional archive-zip SHA-256 pin (enforced if set)
DATASET_TREE_SHA256 = None                       # optional extracted-tree content hash pin (enforced if set)

# Dynamic, per-run acquisition provenance, populated by download_dataset() to
# record exactly HOW the dataset used by this run was obtained and what was (and
# was not) verified. Consumed by fnirs_kappa_realdata_v2.py so the frozen summary
# tells the truth about the specific run rather than the best-case download path.
ACQUISITION_PROVENANCE = None


def _sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def dataset_tree_sha256(root: Path) -> str:
    """Deterministic content hash of an EXTRACTED dataset directory.

    Hashes the sorted list of (POSIX relative path, per-file SHA-256) over every
    regular file under ``root``, so the digest is independent of filesystem order,
    mtimes, and the absolute location. This gives a reproducible fingerprint of the
    dataset content even when it is supplied as an already-extracted directory
    (where no archive-zip checksum is available), addressing the case where a
    pre-extracted tree would otherwise be accepted without any content verification.
    """
    root = Path(root)
    entries = []
    for p in sorted(root.rglob("*")):
        if p.is_file():
            rel = p.relative_to(root).as_posix()
            entries.append(f"{rel}:{_sha256_file(p)}")
    blob = "\n".join(entries).encode()
    return hashlib.sha256(blob).hexdigest()


def _md5_file(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _resolve_zenodo_archive():
    """Return (url, md5) of the first .zip file in the pinned Zenodo record, or
    None if the API cannot be reached."""
    import json as _json
    api = f"https://zenodo.org/api/records/{DATASET_ZENODO_RECORD}"
    try:
        with urllib.request.urlopen(api, timeout=60) as r:
            rec = _json.load(r)
        for f in rec.get("files", []):
            key = f.get("key", "")
            if key.endswith(".zip"):
                url = f.get("links", {}).get("self") or f.get("links", {}).get("download")
                md5 = (f.get("checksum", "") or "").replace("md5:", "")
                if url:
                    return url, md5
    except Exception as e:
        print(f"  [warning] could not reach Zenodo API ({e})")
    return None


def download_dataset(dest: Path) -> Path:
    """Download the BIDS-NIRS-Tapping dataset from the pinned, immutable Zenodo
    release (DATASET_DOI) if not already present, verifying the Zenodo-published
    md5 checksum so the input data is reproducible.

    Parameters
    ----------
    dest : Path
        Target directory for the extracted dataset.

    Returns
    -------
    Path
        Path to the dataset root (containing ``sub-*`` directories).
    """
    global ACQUISITION_PROVENANCE
    if dest.exists() and any(dest.glob("sub-*")):
        # An already-extracted tree carries no archive-zip checksum, so verify a
        # deterministic CONTENT hash instead of accepting it blindly. Enforce the
        # pin if one is set; otherwise record/print it so it can be frozen once.
        # NOTE: this path does NOT contact Zenodo, so no Zenodo MD5 check is
        # performed here; the tree hash fingerprints exactly what was analyzed but
        # is an authenticity check only when compared against a pinned digest.
        tree = dataset_tree_sha256(dest)
        pin_verified = False
        if DATASET_TREE_SHA256 is not None:
            if tree != DATASET_TREE_SHA256:
                raise ValueError(
                    f"extracted-dataset content hash mismatch (expected "
                    f"{DATASET_TREE_SHA256[:12]}..., got {tree[:12]}...). Refusing to "
                    f"proceed with an unverified pre-extracted dataset.")
            pin_verified = True
        ACQUISITION_PROVENANCE = dict(
            acquisition_method="pre-extracted directory (no Zenodo download this run)",
            zenodo_record=DATASET_ZENODO_RECORD, zenodo_doi=DATASET_DOI,
            zenodo_md5_verified=False,
            archive_sha256_verified=False,
            content_tree_sha256=tree,
            tree_hash_verified_against_pin=pin_verified)
        print(f"  Dataset already present at {dest} (no Zenodo download this run)")
        print(f"  dataset content (tree) SHA-256 = {tree}"
              + ("" if DATASET_TREE_SHA256 is None else " (verified against pin)"))
        return dest

    resolved = _resolve_zenodo_archive()
    if resolved is None:
        raise RuntimeError(
            f"Could not resolve the pinned BIDS-NIRS-Tapping archive from Zenodo "
            f"(DOI {DATASET_DOI}). Download it manually from "
            f"https://doi.org/{DATASET_DOI}, extract it so that '{dest}' contains the "
            f"sub-* directories, and re-run. (We deliberately do not fall back to the "
            f"mutable GitHub master branch.)")
    url, md5 = resolved
    zip_path = dest.parent / "bids_nirs_tapping.zip"
    print(f"  Downloading pinned dataset (Zenodo {DATASET_DOI}) from {url} ...")
    urllib.request.urlretrieve(url, zip_path)

    md5_verified = False
    if md5:
        got = _md5_file(zip_path)
        if got != md5:
            zip_path.unlink(missing_ok=True)
            raise ValueError(
                f"BIDS-NIRS-Tapping archive md5 mismatch: Zenodo lists {md5[:12]}..., "
                f"got {got[:12]}.... Refusing to proceed with an unverified input.")
        md5_verified = True
        print(f"  Zenodo md5 verified: {md5[:12]}...")
    archive_sha_verified = False
    if DATASET_SHA256 is not None:
        sha = _sha256_file(zip_path)
        if sha != DATASET_SHA256:
            zip_path.unlink(missing_ok=True)
            raise ValueError(f"archive SHA-256 mismatch (expected {DATASET_SHA256[:12]}..., "
                             f"got {sha[:12]}...).")
        archive_sha_verified = True
    else:
        print(f"  archive SHA-256 = {_sha256_file(zip_path)} (optionally pin via DATASET_SHA256)")

    print("  Extracting ...")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(dest.parent)

    # Zenodo archives of a GitHub repo extract to a top-level dir; find it.
    cands = [p for p in dest.parent.glob("*BIDS-NIRS-Tapping*") if p.is_dir() and p != dest]
    if cands:
        cands[0].rename(dest)
    zip_path.unlink(missing_ok=True)
    tree = dataset_tree_sha256(dest)
    pin_verified = bool(DATASET_TREE_SHA256 is not None and tree == DATASET_TREE_SHA256)
    if DATASET_TREE_SHA256 is not None and not pin_verified:
        raise ValueError(f"extracted-tree hash mismatch after download (expected "
                         f"{DATASET_TREE_SHA256[:12]}..., got {tree[:12]}...).")
    ACQUISITION_PROVENANCE = dict(
        acquisition_method="downloaded pinned Zenodo archive",
        zenodo_record=DATASET_ZENODO_RECORD, zenodo_doi=DATASET_DOI,
        zenodo_md5_verified=md5_verified,
        archive_sha256_verified=archive_sha_verified,
        content_tree_sha256=tree,
        tree_hash_verified_against_pin=pin_verified)
    print(f"  Dataset ready at {dest} (pinned Zenodo {DATASET_DOI})")
    return dest


# ========================================================================
#  2. Data loading
# ========================================================================

def load_subject(bids_root: Path, subject: str) -> mne.io.Raw:
    """Load one subject from the BIDS-NIRS dataset.

    Parameters
    ----------
    bids_root : Path
        Root directory of the BIDS dataset.
    subject : str
        Subject identifier (e.g. ``"01"``).

    Returns
    -------
    mne.io.Raw
        Preloaded continuous fNIRS recording.
    """
    bp = BIDSPath(
        subject=subject,
        task="tapping",
        datatype="nirs",
        root=str(bids_root),
        suffix="nirs",
    )
    raw = read_raw_bids(bp, verbose="ERROR")
    raw.load_data()
    return raw


# ========================================================================
#  3. Channel classification
# ========================================================================

def split_channels(
    raw: mne.io.Raw,
) -> tuple[dict[int, NDArray], dict[int, NDArray], list[int]]:
    """Classify channels into long and short separation by wavelength.

    Parameters
    ----------
    raw : mne.io.Raw
        fNIRS recording whose channel names follow the ``"S1_D1 760"``
        convention.

    Returns
    -------
    long_picks : dict[int, NDArray]
        Channel indices with SDS >= ``SS_THRESHOLD``, keyed by wavelength.
    short_picks : dict[int, NDArray]
        Channel indices with SDS < ``SS_THRESHOLD``, keyed by wavelength.
    wavelengths : list[int]
        Sorted list of wavelengths present in the data.
    """
    dists = source_detector_distances(raw.info)

    long_picks: dict[int, list[int]] = {}
    short_picks: dict[int, list[int]] = {}
    wavelengths_found: set[int] = set()

    for idx, name in enumerate(raw.ch_names):
        parts = name.split()
        if len(parts) < 2:
            continue
        try:
            wl = int(float(parts[-1]))
        except ValueError:
            continue

        wavelengths_found.add(wl)
        bucket = short_picks if dists[idx] < SS_THRESHOLD else long_picks
        bucket.setdefault(wl, []).append(idx)

    # Convert lists to arrays
    long_out = {k: np.array(v) for k, v in long_picks.items()}
    short_out = {k: np.array(v) for k, v in short_picks.items()}

    return long_out, short_out, sorted(wavelengths_found)


# ========================================================================
#  4. Signal processing
# ========================================================================

def bandpass_filter(
    data: NDArray,
    fs: float,
    lo: float = 0.01,
    hi: float = 0.1,
    order: int = 4,
) -> NDArray:
    """Apply a zero-phase Butterworth bandpass filter.

    Parameters
    ----------
    data : NDArray, shape (n_samples, n_channels)
        Time-series data.
    fs : float
        Sampling frequency in Hz.
    lo, hi : float
        Lower and upper cut-off frequencies in Hz.
    order : int
        Filter order.

    Returns
    -------
    NDArray
        Filtered data, same shape as input.
    """
    sos = sp_signal.butter(order, [lo, hi], btype="band", fs=fs, output="sos")
    return sp_signal.sosfiltfilt(sos, data, axis=0)


def ssr_per_wavelength(
    od_long: NDArray, od_short: NDArray
) -> tuple[NDArray, NDArray, NDArray]:
    """Per-wavelength short-separation regression at the OD level.

    For each long channel, the mean short-channel OD signal is regressed
    out via ordinary least squares.

    Parameters
    ----------
    od_long : NDArray, shape (T, n_long)
        Optical-density time series for long-separation channels at a
        single wavelength.
    od_short : NDArray, shape (T, n_short)
        Optical-density time series for short-separation channels at the
        same wavelength.

    Returns
    -------
    od_corrected : NDArray, shape (T, n_long)
        SSR-corrected long-channel OD.
    betas : NDArray, shape (n_long,)
        Regression coefficients.
    r_squared : NDArray, shape (n_long,)
        Coefficient of determination (R^2) for each channel, representing
        the fraction of long-channel variance explained by the short-
        separation regressor.
    """
    ss_mean = od_short.mean(axis=1)                  # (T,)
    ss_sumsq = float(ss_mean @ ss_mean)              # sum of squares

    if ss_sumsq < 1e-20:
        return (od_long.copy(), np.zeros(od_long.shape[1]),
                np.zeros(od_long.shape[1]))

    n_long = od_long.shape[1]
    betas = np.zeros(n_long)
    r_squared = np.zeros(n_long)
    corrected = np.empty_like(od_long)

    for ch in range(n_long):
        beta = float(ss_mean @ od_long[:, ch]) / ss_sumsq
        residual = od_long[:, ch] - beta * ss_mean
        corrected[:, ch] = residual
        betas[ch] = beta

        # R^2 = 1 - Var(residual) / Var(original).
        # Note: uses variance (centred) rather than raw sum-of-squares.
        # For bandpass-filtered signals (zero mean) the two definitions
        # are numerically equivalent; empirically the difference is <0.1 %.
        var_orig = float(np.var(od_long[:, ch]))
        var_resid = float(np.var(residual))
        r_squared[ch] = 1.0 - var_resid / var_orig if var_orig > 1e-20 else 0.0

    return corrected, betas, r_squared


def compute_kappa_ssr(r_squared: NDArray) -> NDArray:
    """Compute the per-channel V_SSR diagnostic from R^2.

    V_SSR = 1 / (1 - R^2_SS) is an inverse residual-variance ratio reported for
    diagnostic purposes only; it is NOT an amplitude-restoration factor and is
    not applied to the data. R^2_SS is the fraction of long-channel variance
    removed by (coupled to) the short-channel regression.

    Parameters
    ----------
    r_squared : NDArray, shape (n_channels,)
        Fraction of long-channel variance explained by SSR.

    Returns
    -------
    v_ssr : NDArray, shape (n_channels,)
        Per-channel V_SSR diagnostic values. Clamped to [1, 50] to
        avoid numerical instability when R^2 approaches 1. Reported only,
        never applied.
    """
    kappa = 1.0 / (1.0 - np.clip(r_squared, 0.0, 0.98))
    return np.clip(kappa, 1.0, 50.0)


# ========================================================================
#  5. MBLL inversion
# ========================================================================

def mbll_inversion(
    od_wl1: NDArray,
    od_wl2: NDArray,
    wl1: int,
    wl2: int,
    sds_mm: float,
) -> tuple[NDArray, NDArray]:
    r"""Solve the two-wavelength MBLL system for HbO2 and HbR.

    The system of equations is:

    .. math::

        \Delta\mathrm{OD}(\lambda) =
            \bigl[\varepsilon_{\mathrm{HbO_2}}(\lambda)\,\Delta c_{\mathrm{HbO_2}}
            + \varepsilon_{\mathrm{HbR}}(\lambda)\,\Delta c_{\mathrm{HbR}}\bigr]
            \cdot d \cdot \mathrm{DPF}(\lambda)

    where *d* is the source-detector separation. Extinction coefficients
    are converted from cm-1 / mM to mm-1 / uM so that the output
    concentrations are in micromolar.

    Parameters
    ----------
    od_wl1, od_wl2 : NDArray, shape (T,)
        Optical-density time series at wavelengths *wl1* and *wl2*.
    wl1, wl2 : int
        Wavelengths in nanometres.
    sds_mm : float
        Source-detector separation in millimetres.

    Returns
    -------
    hbo2, hbr : NDArray, shape (T,)
        Concentration changes in micromolar.
    """
    # Convert MNE's natural-log OD (ln I/I0) to decadic absorbance so it is
    # consistent with the decadic (base-10) extinction coefficients below.
    od_wl1 = od_wl1 / LN10
    od_wl2 = od_wl2 / LN10

    # Normalise OD by optical pathlength
    b1 = od_wl1 / (sds_mm * DPF[wl1])
    b2 = od_wl2 / (sds_mm * DPF[wl2])

    # Extinction matrix in mm-1 / uM
    E = np.array([
        [EXTINCTION_COEFFICIENTS[wl1]["HbO2"] * UNIT_CONVERSION,
         EXTINCTION_COEFFICIENTS[wl1]["HbR"] * UNIT_CONVERSION],
        [EXTINCTION_COEFFICIENTS[wl2]["HbO2"] * UNIT_CONVERSION,
         EXTINCTION_COEFFICIENTS[wl2]["HbR"] * UNIT_CONVERSION],
    ])
    Einv = np.linalg.inv(E)

    hbo2 = Einv[0, 0] * b1 + Einv[0, 1] * b2
    hbr = Einv[1, 0] * b1 + Einv[1, 1] * b2
    return hbo2, hbr


# ========================================================================
#  6. Event handling and epoching
# ========================================================================

def get_tapping_events(raw: mne.io.Raw) -> NDArray:
    """Extract tapping-onset events from the recording's annotations.

    The function identifies tapping conditions by annotation label
    (matching ``"tapping"``, ``"tap"``, or common numeric codes).

    Parameters
    ----------
    raw : mne.io.Raw
        Annotated fNIRS recording.

    Returns
    -------
    NDArray, shape (n_events, 3)
        MNE-format event array (sample, 0, event_id).
    """
    events, event_id = mne.events_from_annotations(raw, verbose="ERROR")

    tapping_ids: list[int] = []
    for name, eid in event_id.items():
        name_lower = name.lower()
        if "tapping" in name_lower or "tap" in name_lower:
            tapping_ids.append(eid)
        if name_lower in ("1.0", "2.0", "1", "2"):
            tapping_ids.append(eid)

    if not tapping_ids:
        all_ids = sorted(event_id.values())
        tapping_ids = all_ids[:-1] if len(all_ids) > 1 else all_ids

    mask = np.isin(events[:, 2], tapping_ids)
    tapping_events = events[mask]

    print(f"  Found {len(tapping_events)} tapping events "
          f"(event IDs: {sorted(set(tapping_events[:, 2]))})")

    return tapping_events


def block_average(
    signal_1d: NDArray,
    event_samples: NDArray,
    fs: float,
    pre: float = 5.0,
    post: float = 20.0,
) -> tuple[NDArray, NDArray, NDArray, int]:
    """Compute a baseline-corrected block average of a 1-D signal.

    Parameters
    ----------
    signal_1d : NDArray, shape (T,)
        Continuous time series.
    event_samples : NDArray
        Onset indices (in samples) for each block.
    fs : float
        Sampling frequency in Hz.
    pre, post : float
        Seconds before and after onset to include.

    Returns
    -------
    time : NDArray
        Time axis relative to onset (seconds).
    mean : NDArray
        Mean across blocks.
    std : NDArray
        Standard deviation across blocks.
    n_good : int
        Number of blocks included.
    """
    pre_samp = int(pre * fs)
    post_samp = int(post * fs)
    n = pre_samp + post_samp

    epochs: list[NDArray] = []
    for s in event_samples:
        lo = s - pre_samp
        hi = s + post_samp
        if lo < 0 or hi > len(signal_1d):
            continue
        epoch = signal_1d[lo:hi].copy()
        epoch -= epoch[:pre_samp].mean()  # baseline correction
        epochs.append(epoch)

    if not epochs:
        return np.zeros(n), np.zeros(n), np.zeros(n), 0

    stacked = np.array(epochs)
    time = np.linspace(-pre, post - 1.0 / fs, n)
    return time, stacked.mean(axis=0), stacked.std(axis=0), len(epochs)


# Physiological response window (seconds after task onset) used for peak
# extraction.  The block-averaged hemodynamic response is searched for its
# peak magnitude only within this window, rather than over the whole
# -pre..post epoch, so that a stray excursion in the pre-stimulus baseline or
# the late post-stimulus tail cannot be selected as the "peak".  The window
# 2-15 s brackets the canonical HbO2 peak (~5-8 s) and the later HbR extremum
# (~8-12 s) for the 5 s tapping block used here.
PEAK_WINDOW_S: tuple[float, float] = (2.0, 15.0)


def windowed_peak_abs(
    mean: NDArray,
    time: NDArray,
    window: tuple[float, float] = PEAK_WINDOW_S,
) -> float:
    """Peak absolute amplitude of a block-averaged response within a window.

    Restricts the max-|amplitude| search to ``window`` (seconds relative to
    onset).  Falls back to the full trace if the window selects no samples.
    """
    t0, t1 = window
    sel = (time >= t0) & (time <= t1)
    seg = mean[sel] if np.any(sel) else mean
    return float(np.max(np.abs(seg)))


# ========================================================================
#  7. Visualisation
# ========================================================================

def generate_figures(
    results: dict[str, Any],
    save_dir: str,
    subject_id: str,
    f_cortex: dict[int, float],
) -> None:
    """Generate and save three publication-quality figures.

    Figures produced:

    1. **Time series** -- Full recording of uncorrected vs. corrected
       HbO2 and HbR, with tapping blocks shaded.
    2. **Hemodynamic response function (HRF)** -- Block-averaged response
       with mean +/- SD envelopes.
    3. **Summary panel** -- Bar chart of peak amplitudes and a text
       summary of correction parameters.

    Parameters
    ----------
    results : dict
        Output dictionary from the main pipeline.
    save_dir : str
        Directory in which to save PNGs.
    subject_id : str
        Subject label for figure titles.
    f_cortex : dict[int, float]
        Cortical sensitivity fractions used in correction.
    """
    plt.rcParams.update(_PLOT_PARAMS)

    t = results["time"]
    ba = results["block_avg"]

    # ── Figure 1: full time series ──────────────────────────────────
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)
    fig.suptitle(
        f"Time Series — Uncorrected vs κ-Corrected fNIRS\n"
        f"(Subject {subject_id}, SDS ≈ {results['sds_mm']:.0f} mm)",
        fontsize=13, fontweight="bold",
    )

    ax = axes[0]
    ax.plot(t, results["hbo2_uncorr"], "--", color="#FF6B6B", lw=1.2,
            alpha=0.7, label="Uncorrected MBLL")
    ax.plot(t, results["hbo2_corr"], "-", color="#C92A2A", lw=1.8,
            label="κ-Corrected")
    ax.set_ylabel("Δ[HbO₂] (μM)")
    ax.legend(loc="upper right", frameon=False)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    ax.plot(t, results["hbr_uncorr"], "--", color="#74C0FC", lw=1.2,
            alpha=0.7, label="Uncorrected MBLL")
    ax.plot(t, results["hbr_corr"], "-", color="#1971C2", lw=1.8,
            label="κ-Corrected")
    ax.set_ylabel("Δ[HbR] (μM)")
    ax.set_xlabel("Time (s)")
    ax.legend(loc="upper right", frameon=False)
    ax.grid(True, alpha=0.3)

    for onset_s in results["event_times"]:
        for a in axes:
            a.axvspan(onset_s, onset_s + 5.0, color="gray",
                      alpha=0.08, zorder=0)

    plt.tight_layout()
    path = os.path.join(save_dir, "Figure_RealData_TimeSeries.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")

    # ── Figure 2: block-averaged HRF ────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Block-Averaged Hemodynamic Response — Real fNIRS Data\n"
        f"(mean ± SD across {ba['n_blocks']} tapping blocks)",
        fontsize=13, fontweight="bold",
    )

    for ax, chrom, unit_label in zip(
        axes,
        [("hbo2_uncorr", "hbo2_corr"), ("hbr_uncorr", "hbr_corr")],
        ["Δ[HbO₂] (μM)", "Δ[HbR] (μM)"],
    ):
        t_ba = ba["time"]
        m_unc = ba[chrom[0] + "_mean"]
        s_unc = ba[chrom[0] + "_std"]
        m_cor = ba[chrom[1] + "_mean"]
        s_cor = ba[chrom[1] + "_std"]

        c_unc = "#FF6B6B" if "hbo2" in chrom[0] else "#74C0FC"
        c_cor = "#C92A2A" if "hbo2" in chrom[0] else "#1971C2"

        ax.fill_between(t_ba, m_unc - s_unc, m_unc + s_unc,
                        color=c_unc, alpha=0.15)
        ax.plot(t_ba, m_unc, "--", color=c_unc, lw=1.5,
                label="Uncorrected")
        ax.fill_between(t_ba, m_cor - s_cor, m_cor + s_cor,
                        color=c_cor, alpha=0.15)
        ax.plot(t_ba, m_cor, "-", color=c_cor, lw=2.0,
                label="κ-Corrected")

        ax.axvspan(0, 5.0, color="gray", alpha=0.10, label="Task")
        ax.set_xlabel("Time relative to onset (s)")
        ax.set_ylabel(unit_label)
        ax.legend(frameon=False)
        ax.grid(True, alpha=0.3)

    plt.tight_layout()
    path = os.path.join(save_dir, "Figure_RealData_HRF.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")

    # ── Figure 3: summary panel ─────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle(
        "κ-Correction Summary — Real fNIRS Data\n"
        f"(Subject {subject_id}, SDS ≈ {results['sds_mm']:.0f} mm)",
        fontsize=13, fontweight="bold",
    )

    hbo2_peak_unc = windowed_peak_abs(ba["hbo2_uncorr_mean"], ba["time"])
    hbo2_peak_cor = windowed_peak_abs(ba["hbo2_corr_mean"], ba["time"])
    hbr_peak_unc = windowed_peak_abs(ba["hbr_uncorr_mean"], ba["time"])
    hbr_peak_cor = windowed_peak_abs(ba["hbr_corr_mean"], ba["time"])

    x = np.arange(2)
    w = 0.35
    ax = axes[0]
    ax.bar(x - w / 2, [hbo2_peak_unc, hbo2_peak_cor], w,
           color=["#FFB3B3", "#C92A2A"], label="HbO₂")
    ax.bar(x + w / 2, [hbr_peak_unc, hbr_peak_cor], w,
           color=["#A9DDFF", "#1971C2"], label="HbR")
    ax.set_xticks(x)
    ax.set_xticklabels(["Uncorrected", "κ-Corrected"])
    ax.set_ylabel("Peak |Amplitude| (μM)")
    ax.set_title("Block-Averaged Peak Concentrations")
    ax.legend(frameon=False)
    ax.grid(True, axis="y", alpha=0.3)

    # Net scaling ratios
    ratio_hbo2 = hbo2_peak_cor / hbo2_peak_unc if hbo2_peak_unc > 0 else 0
    ratio_hbr = hbr_peak_cor / hbr_peak_unc if hbr_peak_unc > 0 else 0

    wl1, wl2 = results["wavelengths"]
    k_ssr_wl = results.get("kappa_ssr_by_wl", {})
    summary_text = (
        f"Subject: {subject_id}\n"
        f"SDS (long): {results['sds_mm']:.1f} mm\n"
        f"Wavelengths: {wl1}, {wl2} nm\n"
        f"Sampling rate: {results['fs']:.1f} Hz\n"
        f"Tapping blocks: {ba['n_blocks']}\n"
        f"\n"
        f"f_cortex({wl1}): {f_cortex.get(wl1, 'N/A')}\n"
        f"f_cortex({wl2}): {f_cortex.get(wl2, 'N/A')}\n"
        f"κ_PV (mean, applied):    {results['kappa_pv_mean']:.3f}\n"
        f"V_SSR({wl1}) [diag.]:    {k_ssr_wl.get(wl1, 'N/A'):.3f}\n"
        f"V_SSR({wl2}) [diag.]:    {k_ssr_wl.get(wl2, 'N/A'):.3f}\n"
        f"\n"
        f"Peak |HbO₂| uncorr: {hbo2_peak_unc:.4f} μM\n"
        f"Peak |HbO₂| corr:   {hbo2_peak_cor:.4f} μM\n"
        f"Net scaling:         {ratio_hbo2:.2f}×\n"
        f"\n"
        f"Peak |HbR| uncorr:   {hbr_peak_unc:.4f} μM\n"
        f"Peak |HbR| corr:     {hbr_peak_cor:.4f} μM\n"
        f"Net scaling:         {ratio_hbr:.2f}×\n"
    )

    ax = axes[1]
    ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
            fontsize=10, va="top", family="monospace",
            bbox=dict(boxstyle="round", fc="wheat", alpha=0.3))
    ax.axis("off")
    ax.set_title("Correction Summary")

    plt.tight_layout()
    path = os.path.join(save_dir, "Figure_RealData_Summary.png")
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {path}")


# ========================================================================
#  Main pipeline
# ========================================================================

def main() -> None:
    """Execute the full kappa-correction pipeline on real fNIRS data."""

    print("=" * 70)
    print("  In Vivo Proof-of-Concept: κ-Correction of fNIRS Data")
    print("=" * 70)

    script_dir = Path(__file__).resolve().parent
    data_dir = script_dir / "BIDS-NIRS-Tapping-data"

    # ── Step 1: dataset acquisition ─────────────────────────────────
    print("\n[1] Dataset acquisition")
    data_dir = download_dataset(data_dir)

    # ── Step 2: load subject ────────────────────────────────────────
    subject = "01"
    print(f"\n[2] Loading subject {subject} ...")
    try:
        raw = load_subject(data_dir, subject)
    except Exception as exc:
        print(f"  ERROR loading subject: {exc}")
        print("\n  Troubleshooting:")
        print('    pip install "numpy<2" h5py --force-reinstall')
        sys.exit(1)

    fs = raw.info["sfreq"]
    print(f"  Channels: {len(raw.ch_names)}")
    print(f"  Sampling rate: {fs} Hz")
    print(f"  Duration: {raw.times[-1]:.1f} s")

    # ── Step 3: channel classification ──────────────────────────────
    print("\n[3] Classifying channels by wavelength and separation ...")
    long_picks, short_picks, wavelengths = split_channels(raw)
    print(f"  Wavelengths: {wavelengths}")
    for wl in wavelengths:
        n_long = len(long_picks.get(wl, []))
        n_short = len(short_picks.get(wl, []))
        print(f"    {wl} nm -- {n_long} long, {n_short} short")

    if len(wavelengths) < 2:
        print("ERROR: At least two wavelengths are required.")
        sys.exit(1)

    wl1, wl2 = wavelengths[0], wavelengths[1]

    # Ensure physical constants exist for the data wavelengths
    for wl in (wl1, wl2):
        if wl not in EXTINCTION_COEFFICIENTS:
            nearest = min(EXTINCTION_COEFFICIENTS, key=lambda x: abs(x - wl))
            EXTINCTION_COEFFICIENTS[wl] = EXTINCTION_COEFFICIENTS[nearest]
        if wl not in DPF:
            nearest = min(DPF, key=lambda x: abs(x - wl))
            DPF[wl] = DPF[nearest]

    # Compute median SDS for long channels
    dists = source_detector_distances(raw.info)
    long_all = np.concatenate([long_picks[wl1], long_picks[wl2]])
    sds_mm = float(np.median(dists[long_all])) * 1000.0
    print(f"  Median long-channel SDS: {sds_mm:.1f} mm")

    # Interpolate cortical sensitivity fractions for this SDS
    f_cortex = get_f_cortex(sds_mm)
    print(f"  f_cortex (interpolated for SDS = {sds_mm:.1f} mm):")
    for wl_key in sorted(f_cortex):
        print(f"    {wl_key} nm: {f_cortex[wl_key]:.4f}")

    # Ensure f_cortex covers both data wavelengths
    for wl in (wl1, wl2):
        if wl not in f_cortex:
            nearest = min(f_cortex, key=lambda x: abs(x - wl))
            f_cortex[wl] = f_cortex[nearest]

    # ── Step 4: intensity to optical density ────────────────────────
    print("\n[4] Converting intensity to optical density ...")
    raw_od = optical_density(raw)
    od_data = raw_od.get_data().T  # (T, n_channels)

    # ── Step 5: bandpass filtering ──────────────────────────────────
    print("[5] Bandpass filtering (0.01 -- 0.1 Hz) ...")
    od_filt = bandpass_filter(od_data, fs, lo=0.01, hi=0.1)

    # ── Step 6: short-separation regression ─────────────────────────
    print("[6] Short-separation regression (per wavelength) ...")
    od_ssr = od_filt.copy()
    r2_by_wl: dict[int, NDArray] = {}
    kappa_ssr_by_wl: dict[int, float] = {}
    for wl in (wl1, wl2):
        lp = long_picks[wl]
        sp = short_picks.get(wl, np.array([], dtype=int))
        if len(sp) == 0:
            print(f"    {wl} nm -- no short channels, skipping SSR")
            r2_by_wl[wl] = np.zeros(len(lp))
            kappa_ssr_by_wl[wl] = 1.0
            continue
        corrected, betas, r2 = ssr_per_wavelength(od_filt[:, lp], od_filt[:, sp])
        od_ssr[:, lp] = corrected
        r2_by_wl[wl] = r2
        k_ssr = compute_kappa_ssr(r2)
        kappa_ssr_by_wl[wl] = float(np.mean(k_ssr))
        print(f"    {wl} nm -- mean beta = {betas.mean():.4f}, "
              f"mean R^2 = {r2.mean():.4f}, "
              f"mean V_SSR = {np.mean(k_ssr):.3f}")

    # ── Step 7: V_SSR is a DIAGNOSTIC, not applied ───────────────
    # V_SSR = 1/(1 - R^2_SS) was previously multiplied into the long-channel
    # OD here.  That double-counts the dilution: short-separation regression has
    # already removed the SUPERFICIAL (systemic) component, and the residual
    # cortical signal is recovered by the partial-volume factor kappa_PV =
    # 1/f_cortex alone (Step 8).  Restoring the SSR-removed variance with
    # V_SSR would re-inject the systemic signal SSR was meant to remove, and
    # together with the corrected (large) kappa_PV drove the corrected amplitude
    # to a non-physiological ~11 uM.  We therefore report V_SSR as a
    # diagnostic of long-channel variance removed by/coupled to the short-channel
    # regression (not of residual superficial contamination) but do NOT apply it,
    # consistent with the synthetic pipeline (where only kappa_PV is applied).
    # See KAPPA_PV_ATLAS_CHECK.md (Rec. 2) and CHANGELOG_FIXES.md.
    print("[7] V_SSR reported as diagnostic (NOT applied; see Step 8) ...")
    od_ssr_corrected = od_ssr.copy()
    for wl in (wl1, wl2):
        k_ssr = compute_kappa_ssr(r2_by_wl[wl])
        print(f"    {wl} nm -- V_SSR (diagnostic) range: "
              f"[{k_ssr.min():.3f}, {k_ssr.max():.3f}], mean {np.mean(k_ssr):.3f}")

    # ── Step 8: partial-volume correction ───────────────────────────
    print("[8] Partial-volume correction ...")
    od_pv = od_ssr_corrected.copy()
    kappa_pvs: dict[int, float] = {}
    for wl in (wl1, wl2):
        f = f_cortex[wl]
        kpv = 1.0 / f
        kappa_pvs[wl] = kpv
        lp = long_picks[wl]
        od_pv[:, lp] = od_ssr_corrected[:, lp] / f
        print(f"    {wl} nm -- f_cortex = {f:.4f}, kappa_PV = {kpv:.3f}")

    kappa_pv_mean = float(np.mean(list(kappa_pvs.values())))
    # No cross-wavelength mean of V_SSR (or kappa_total) is formed: V_SSR
    # is a per-wavelength diagnostic only, and the 760/850 nm values differ too
    # much for their arithmetic mean to be a meaningful quantity.

    # ── Step 9: MBLL inversion ──────────────────────────────────────
    print("[9] MBLL inversion ...")

    def _parse_sd_label(name: str) -> str:
        """Extract source-detector label (e.g. 'S1_D1') from channel name."""
        return name.split()[0] if " " in name else name

    sd_wl1 = {_parse_sd_label(raw.ch_names[i]): i for i in long_picks[wl1]}
    sd_wl2 = {_parse_sd_label(raw.ch_names[i]): i for i in long_picks[wl2]}
    common_labels = sorted(set(sd_wl1) & set(sd_wl2))

    if not common_labels:
        print("ERROR: No matched channel pairs found.")
        sys.exit(1)

    print(f"  Matched channel pairs: {len(common_labels)}")

    n_t = od_filt.shape[0]
    hbo2_uncorr_all: list[NDArray] = []
    hbr_uncorr_all: list[NDArray] = []
    hbo2_corr_all: list[NDArray] = []
    hbr_corr_all: list[NDArray] = []

    for label in common_labels:
        i1, i2 = sd_wl1[label], sd_wl2[label]

        # Uncorrected pathway: filtered OD -> MBLL
        hbo2_u, hbr_u = mbll_inversion(
            od_filt[:, i1], od_filt[:, i2], wl1, wl2, sds_mm
        )
        hbo2_uncorr_all.append(hbo2_u)
        hbr_uncorr_all.append(hbr_u)

        # Corrected pathway: SSR + PV corrected OD -> MBLL
        hbo2_c, hbr_c = mbll_inversion(
            od_pv[:, i1], od_pv[:, i2], wl1, wl2, sds_mm
        )
        hbo2_corr_all.append(hbo2_c)
        hbr_corr_all.append(hbr_c)

    # Channel-averaged traces
    hbo2_uncorr = np.mean(hbo2_uncorr_all, axis=0)
    hbr_uncorr = np.mean(hbr_uncorr_all, axis=0)
    hbo2_corr = np.mean(hbo2_corr_all, axis=0)
    hbr_corr = np.mean(hbr_corr_all, axis=0)

    time = np.arange(n_t) / fs

    # ── Step 10: event extraction and block averaging ───────────────
    print("\n[10] Extracting tapping events and block-averaging ...")
    tapping_events = get_tapping_events(raw)
    event_samples = tapping_events[:, 0]
    event_times = event_samples / fs

    ba: dict[str, Any] = {}
    for name, sig in [
        ("hbo2_uncorr", hbo2_uncorr),
        ("hbr_uncorr", hbr_uncorr),
        ("hbo2_corr", hbo2_corr),
        ("hbr_corr", hbr_corr),
    ]:
        t_ba, m, s, n = block_average(sig, event_samples, fs,
                                       pre=5.0, post=20.0)
        ba["time"] = t_ba
        ba[name + "_mean"] = m
        ba[name + "_std"] = s
        ba["n_blocks"] = n

    print(f"  Blocks used: {ba['n_blocks']}")

    # ── Step 11: results summary ────────────────────────────────────
    # Peaks are taken within the physiological response window (PEAK_WINDOW_S)
    # rather than over the whole epoch, so a baseline/tail excursion cannot be
    # mistaken for the response peak.
    hbo2_peak_unc = windowed_peak_abs(ba["hbo2_uncorr_mean"], ba["time"])
    hbo2_peak_cor = windowed_peak_abs(ba["hbo2_corr_mean"], ba["time"])
    hbr_peak_unc = windowed_peak_abs(ba["hbr_uncorr_mean"], ba["time"])
    hbr_peak_cor = windowed_peak_abs(ba["hbr_corr_mean"], ba["time"])

    print("\n" + "=" * 70)
    print("  RESULTS SUMMARY")
    print("=" * 70)
    print(f"  Subject:              {subject}")
    print(f"  SDS (long):           {sds_mm:.1f} mm")
    print(f"  Wavelengths:          {wl1}, {wl2} nm")
    print(f"  Channel pairs:        {len(common_labels)}")
    print(f"  Tapping blocks:       {ba['n_blocks']}")
    print()
    print("  Correction factors:")
    print(f"    kappa_PV (mean, applied):    {kappa_pv_mean:.3f}")
    print(f"    V_SSR({wl1}) [diagnostic]:    {kappa_ssr_by_wl.get(wl1, 0):.3f}")
    print(f"    V_SSR({wl2}) [diagnostic]:    {kappa_ssr_by_wl.get(wl2, 0):.3f}")
    print("    (V_SSR is per-wavelength and NOT applied; no cross-wavelength")
    print("     mean or kappa_total is reported.)")
    print()
    print(f"  Peak |HbO2| uncorr:   {hbo2_peak_unc:.4f} uM")
    print(f"  Peak |HbO2| corr:     {hbo2_peak_cor:.4f} uM")
    if hbo2_peak_unc > 0:
        print(f"  Net scaling (HbO2):   {hbo2_peak_cor / hbo2_peak_unc:.2f}x")
    print()
    print(f"  Peak |HbR| uncorr:    {hbr_peak_unc:.4f} uM")
    print(f"  Peak |HbR| corr:      {hbr_peak_cor:.4f} uM")
    if hbr_peak_unc > 0:
        print(f"  Net scaling (HbR):    {hbr_peak_cor / hbr_peak_unc:.2f}x")
    print("=" * 70)

    # ── Step 12: figure generation ──────────────────────────────────
    print("\n[11] Generating figures ...")
    results = {
        "time": time,
        "hbo2_uncorr": hbo2_uncorr,
        "hbr_uncorr": hbr_uncorr,
        "hbo2_corr": hbo2_corr,
        "hbr_corr": hbr_corr,
        "event_times": event_times,
        "block_avg": ba,
        "sds_mm": sds_mm,
        "fs": fs,
        "wavelengths": [wl1, wl2],
        "kappa_pv_mean": kappa_pv_mean,
        "kappa_ssr_by_wl": kappa_ssr_by_wl,
    }
    generate_figures(results, str(script_dir), subject, f_cortex)

    print("\nDone.")


if __name__ == "__main__":
    main()
