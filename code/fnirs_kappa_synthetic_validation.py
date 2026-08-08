#!/usr/bin/env python3
"""
fNIRS MBLL Bias Quantification and Two-Layer κ Correction Study

Complete analysis pipeline for evaluating κ correction factors in functional
near-infrared spectroscopy (fNIRS) using a two-layer tissue model.

Pipeline includes:
  - Proper two-layer Kienle SDA fluence (Kienle et al. 1998)
  - Full 3D sensitivity volume integral (x, y, z) with y-symmetry
  - Wavelength-dependent cortical fraction (f_cortex)
  - Wavelength-specific κ(PV) correction at OD level
  - Wavelength-specific forward model in synthetic data generation
  - SSR regression on same-wavelength short-channel OD
  - Per-wavelength V(SSR, λ) = 1/(1 − R²_SS(λ)) reported as a diagnostic only
    (NOT applied; no cross-wavelength mean or κ_total is formed)
  - Exact Student-t p-value via betainc
  - Core κ correction framework (DPF, PV, SSR)
  - HbO₂ and HbR performance evaluation
  - Per-subject variability analysis
  - Robustness to superficial layer thickness
  - Robustness to optical property variation
  - Grid convergence study
  - Computational efficiency benchmarking

Authors: Neth Sagara & Amir Ahsan
Institution: Irvine Valley College
Date: February 2026

Reference:
    Sagara N and Ahsan A (2026). "Partial-volume correction in continuous-wave
    fNIRS: operating regime, a Monte-Carlo convergence caveat, and a reproducible
    pipeline." Neurophotonics (submitted).
"""

# =============================================================================
# IMPORTS
# =============================================================================
from __future__ import annotations
import math
import time
import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import json
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.special import j0, betainc

warnings.filterwarnings('ignore', category=RuntimeWarning)

# NumPy compatibility: np.trapezoid (NumPy>=2.0) vs np.trapz (NumPy<2.0)
_np_trapz = getattr(np, 'trapezoid', None) or np.trapz


def _robustness_from_json(require=None):
    """Load the secondary robustness sweeps from results/robustness_secondary.json
    (produced by mc_robustness_sweeps.py) so Figure 3 and the robustness numbers are
    read DIRECTLY from that versioned artifact rather than from embedded constants.

    The file's payload SHA-256 (provenance.data_sha256) is re-verified on load. In
    RELEASE mode (env FNIRS_RELEASE=1, or require=True) the file is MANDATORY and any
    absence/parse-error/hash-mismatch RAISES, so a release build can never silently
    fall back to the archived constant arrays. Otherwise (development) a failure
    returns None and the caller may use the clearly-labelled archived fallback only if
    the explicit dev flag FNIRS_ALLOW_FALLBACK=1 is set.

    Returns dict(thickness_T, thickness_F, opt, source) or None."""
    import os as _os, hashlib as _hashlib
    if require is None:
        require = _os.environ.get("FNIRS_RELEASE") == "1"
    here = _os.path.dirname(_os.path.abspath(__file__))
    cands = [_os.path.join(here, "..", "results", "robustness_secondary.json"),
             _os.path.join(_os.getcwd(), "results", "robustness_secondary.json"),
             _os.path.join(here, "robustness_secondary.json")]
    found = next((c for c in cands if _os.path.exists(c)), None)
    if found is None:
        if require:
            raise FileNotFoundError(
                "RELEASE mode: results/robustness_secondary.json is required (run "
                "mc_robustness_sweeps.py) and must not be replaced by archived arrays.")
        return None
    try:
        with open(found) as fh:
            d = json.load(fh)
        # re-verify the payload hash (same recipe as mc_robustness_sweeps.py). In
        # RELEASE mode the hash is MANDATORY: a file with provenance.data_sha256
        # removed is rejected rather than silently skipped (matching the production
        # loader, which always requires its payload hash).
        recorded = (d.get("provenance") or {}).get("data_sha256")
        if not recorded:
            if require:
                raise ValueError("RELEASE mode: robustness_secondary.json is missing "
                                 "provenance.data_sha256; refusing to use an unverifiable "
                                 "robustness artifact.")
        else:
            payload = {k: d[k] for k in d if k != "provenance"}
            blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
            got = _hashlib.sha256(blob).hexdigest()
            if got != recorded:
                raise ValueError(f"robustness_secondary.json payload hash mismatch "
                                 f"({got[:12]} != {recorded[:12]})")
        rows = sorted(d["thickness_sweep"]["rows"], key=lambda r: r["sup_thickness_mm"])
        T = [float(r["sup_thickness_mm"]) for r in rows]
        F = [float(r["f_cortex"]) for r in rows]
        osw = d["optical_property_sweep"]
        def _col(key, sub):
            return [float(r[sub]) for r in sorted(osw[key], key=lambda r: r["pct"])]
        pct = [float(r["pct"]) for r in sorted(osw["mua"], key=lambda r: r["pct"])]
        mua_f = _col("mua", "f_cortex"); mua_k = _col("mua", "kappa_pv")
        musp_f = _col("musp", "f_cortex"); musp_k = _col("musp", "kappa_pv")
        i0 = pct.index(0.0)
        opt = {'variations_pct': pct,
               'mua': {'f': mua_f, 'k': mua_k},
               'musp': {'f': musp_f, 'k': musp_k},
               'baseline': {'f': mua_f[i0], 'k': mua_k[i0]}}
        return dict(thickness_T=T, thickness_F=F, opt=opt, source=found)
    except Exception as e:
        if require:
            raise
        print(f"  [warning] could not load/verify robustness_secondary.json ({e})")
        return None


def _synth_provenance(produced_by):
    """Provenance block for synthetic-validation JSON outputs, from the shared
    provenance helper (uniform schema across all package scripts)."""
    import provenance as _prov
    return _prov.provenance(produced_by,
                            input_hashes=_prov.fcortex_production_input(),
                            extra=dict(forward_model=_prov.fcortex_production_input(), seed=42))

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def bandpass_filter(data: np.ndarray, sfreq: float, low_freq: float, high_freq: float) -> np.ndarray:
    """
    Simple FFT-based bandpass filter (replaces scipy.signal.filtfilt).

    Parameters:
        data: Input signal
        sfreq: Sampling frequency (Hz)
        low_freq: Low cutoff frequency (Hz)
        high_freq: High cutoff frequency (Hz)

    Returns:
        Filtered signal
    """
    n = len(data)
    freqs = np.fft.rfftfreq(n, d=1/sfreq)
    fft_data = np.fft.rfft(data)

    # Create bandpass mask with smooth transitions
    mask = np.zeros_like(freqs)
    in_band = (freqs >= low_freq) & (freqs <= high_freq)
    mask[in_band] = 1.0

    transition_width = 0.02  # Hz
    for i, f in enumerate(freqs):
        if low_freq - transition_width < f < low_freq:
            mask[i] = (f - (low_freq - transition_width)) / transition_width
        elif high_freq < f < high_freq + transition_width:
            mask[i] = 1 - (f - high_freq) / transition_width

    filtered_fft = fft_data * mask
    filtered_data = np.fft.irfft(filtered_fft, n=n)
    return filtered_data


def pinv(X: np.ndarray) -> np.ndarray:
    """Pseudoinverse using numpy (replaces scipy.linalg.pinv)."""
    return np.linalg.pinv(X)


def pearsonr(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Pearson correlation coefficient with exact two-tailed Student-t p-value.

    Uses the regularised incomplete beta function to evaluate the exact CDF of
    the t-distribution:

        p = I_{ν/(ν+t²)}(ν/2, 1/2)

    which is the exact two-tailed p-value for H₀: ρ = 0.
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)

    x_mean = np.mean(x)
    y_mean = np.mean(y)
    x_centered = x - x_mean
    y_centered = y - y_mean

    numerator = np.sum(x_centered * y_centered)
    denominator = np.sqrt(np.sum(x_centered**2) * np.sum(y_centered**2))

    if denominator == 0:
        return 0.0, 1.0

    r = float(np.clip(numerator / denominator, -1.0, 1.0))
    n = len(x)
    nu = n - 2  # degrees of freedom

    if abs(r) >= 1.0 or nu <= 0:
        # Perfect correlation or insufficient data
        p = 0.0
    else:
        t_sq = r**2 * nu / (1.0 - r**2)
        # Exact two-tailed p-value via regularised incomplete beta function
        # I_x(a, b) = betainc(a, b, x) in scipy.special
        p = float(betainc(nu / 2.0, 0.5, nu / (nu + t_sq)))

    return r, p


# =============================================================================
# WAVELENGTH-DEPENDENT DPF
# =============================================================================

def scholkmann_dpf(wavelength_nm: float, age_years: float = 25.0) -> float:
    """
    Compute the Differential Pathlength Factor using the general equation
    from Scholkmann and Wolf (2013).

    DPF(λ, A) = 223.3 + 0.05624·A^0.8493 − 5.723e-7·λ³ + 0.001245·λ² − 0.9025·λ

    Parameters:
        wavelength_nm: Wavelength in nm (valid range ~690–900 nm)
        age_years: Subject age in years (default 25)

    Returns:
        DPF value (dimensionless)

    Reference:
        Scholkmann F and Wolf M 2013 J. Biomed. Opt. 18 105004
    """
    lam = float(wavelength_nm)
    A = float(age_years)
    dpf = (223.3
           + 0.05624 * A**0.8493
           - 5.723e-7 * lam**3
           + 0.001245 * lam**2
           - 0.9025 * lam)
    return max(dpf, 1.0)  # Guard against non-physical values


def get_dpf_dict(wavelengths: List[int], age_years: float = 25.0) -> Dict[int, float]:
    """
    Compute wavelength-specific DPF values for a list of wavelengths.

    Parameters:
        wavelengths: List of wavelengths (nm)
        age_years: Subject age in years

    Returns:
        Dictionary mapping wavelength → DPF
    """
    return {wl: scholkmann_dpf(wl, age_years) for wl in wavelengths}


# =============================================================================
# DATA CLASSES
# =============================================================================

@dataclass
class ExtinctionCoefficients:
    """
    Decadic (base-10) molar extinction coefficients for HbO2 and HbR
    at common fNIRS wavelengths.

    Units: (cm^-1) / (mM), decadic
    Source: Prahl/Oregon Medical Laser Center (Prahl 1999)
    """
    coefficients: Dict[int, Tuple[float, float]] = field(default_factory=lambda: {
        690: (0.276, 2.052),
        700: (0.290, 1.794),
        730: (0.390, 1.102),
        750: (0.518, 1.405),
        760: (0.586, 1.549),
        780: (0.710, 1.075),
        800: (0.816, 0.762),
        830: (0.974, 0.693),
        850: (1.058, 0.691),
        880: (1.154, 0.726),
    })

    def get_matrix(self, wavelengths: List[int]) -> np.ndarray:
        """Build extinction coefficient matrix E for given wavelengths."""
        E = np.array([self.coefficients[w] for w in wavelengths])
        return E


@dataclass
class OpticalProperties:
    """
    Optical properties for a tissue layer at a specific wavelength.

    Attributes:
        mua: Absorption coefficient (1/mm)
        musp: Reduced scattering coefficient (1/mm)
        n: Refractive index (dimensionless)
    """
    mua: float
    musp: float
    n: float = 1.4

    @property
    def D(self) -> float:
        """Diffusion coefficient (mm)"""
        return 1.0 / (3.0 * (self.mua + self.musp))

    @property
    def mu_eff(self) -> float:
        """Effective attenuation coefficient (1/mm)"""
        return math.sqrt(3.0 * self.mua * (self.mua + self.musp))


@dataclass
class TwoLayerModel:
    """
    Two-layer slab model for head tissue.

    Layer 1 (Superficial): Scalp + Skull + CSF (~10-15mm thick)
    Layer 2 (Cortex): Gray matter (semi-infinite)
    """
    z_superficial_mm: float
    props_superficial: Dict[int, OpticalProperties] = field(default_factory=dict)
    props_cortex: Dict[int, OpticalProperties] = field(default_factory=dict)

    @classmethod
    def default_adult(cls, z_superficial_mm: float = 12.0) -> 'TwoLayerModel':
        """Create a two-layer model with typical adult optical properties."""
        model = cls(z_superficial_mm=z_superficial_mm)
        wavelengths = [690, 730, 760, 780, 800, 830, 850]

        for wl in wavelengths:
            mua_sup = 0.012 + 0.002 * (wl - 700) / 150
            musp_sup = 1.2 - 0.3 * (wl - 700) / 150
            model.props_superficial[wl] = OpticalProperties(mua=mua_sup, musp=musp_sup)

            mua_ctx = 0.018 + 0.003 * (wl - 700) / 150
            musp_ctx = 0.9 - 0.2 * (wl - 700) / 150
            model.props_cortex[wl] = OpticalProperties(mua=mua_ctx, musp=musp_ctx)

        return model


# =============================================================================
# LIGHT TRANSPORT - TWO-LAYER KIENLE SDA FLUENCE
# =============================================================================

def fluence_two_layer_kienle(
    rho_mm: float,
    z_mm: float,
    z_sup_mm: float,
    props_sup: OpticalProperties,
    props_cortex: OpticalProperties,
    n_s: int = 500,
    s_max: float = 30.0
) -> float:
    """
    Compute CW diffusion fluence for two-layer medium using Kienle et al. 1998
    steady-state diffusion approximation via Hankel transform.

    The fluence is obtained via zeroth-order Hankel transform:
        Φ(ρ, z) = (1/2π) ∫₀^∞ phi_k(s, z) · J₀(s·ρ) · s ds

    Parameters:
        rho_mm: Radial distance from source axis (mm)
        z_mm: Depth below surface (mm)
        z_sup_mm: Thickness of superficial layer (mm)
        props_sup: Optical properties of superficial layer
        props_cortex: Optical properties of cortex layer
        n_s: Number of Hankel quadrature points (default 500)
        s_max: Maximum spatial frequency (mm^-1, default 30)

    Returns:
        Fluence (arbitrary units)

    Reference:
        Kienle A, Patterson MS, Dögnitz N, Bays R, Wagnières G and
        van den Bergh H 1998. Noninvasive determination of the optical
        properties of two-layered turbid media. Appl. Opt. 37: 779-91.
    """
    # Quadrature setup
    s_vals = np.linspace(0.001, s_max, n_s)
    ds = (s_max - 0.001) / (n_s - 1)

    # Extract optical properties
    D1 = props_sup.D
    D2 = props_cortex.D
    mua1 = props_sup.mua
    mua2 = props_cortex.mua
    musp1 = props_sup.musp

    # Extrapolated boundary (A ≈ 2.95 for n=1.4)
    A = 2.95
    zb = 2.0 * A * D1

    # Source depth (isotropic)
    zp = 1.0 / musp1

    # Compute Hankel domain solution for all s at once (vectorized)
    alpha1 = np.sqrt((D1 * s_vals**2 + mua1) / D1)
    alpha2 = np.sqrt((D2 * s_vals**2 + mua2) / D2)
    Da = (D1*alpha1 - D2*alpha2) / (D1*alpha1 + D2*alpha2)

    if z_mm < zp:
        arg = zp - z_mm
    else:
        arg = z_mm - zp

    dum1 = np.exp(-alpha1*arg) - np.exp(-alpha1*(2*zb + zp + z_mm))
    dum2 = (np.exp(-alpha1*(-zp + 2*z_sup_mm - z_mm))
            - np.exp(-alpha1*(2*zb + zp + 2*z_sup_mm - z_mm))
            - np.exp(-alpha1*(-zp + 2*z_sup_mm + 2*zb + z_mm))
            + np.exp(-alpha1*(4*zb + zp + 2*z_sup_mm + z_mm)))

    if z_mm < z_sup_mm:
        # Top layer (z < d_sup)
        phi_k = (dum1 + Da*dum2) / (2*D1*alpha1)
    else:
        # Bottom layer (z >= d_sup)
        dum3 = np.exp(alpha2*(z_sup_mm - z_mm)) / (D1*alpha1 + D2*alpha2)
        dum4 = np.exp(alpha1*(zp - z_sup_mm)) - np.exp(-alpha1*(2*zb + zp + z_sup_mm))
        dum5 = (np.exp(alpha1*(zp - 3*z_sup_mm - 2*zb))
                - np.exp(-alpha1*(4*zb + zp + 3*z_sup_mm)))
        phi_k = dum3 * (dum4 - Da*dum5)

    # Inverse Hankel transform to real space
    j0_vals = j0(s_vals * rho_mm)
    integrand = phi_k * j0_vals * s_vals
    fluence = (1.0 / (2.0 * np.pi)) * _np_trapz(integrand, dx=ds)

    return max(float(fluence), 1e-20)


def fluence_two_layer_kienle_vectorized(
    rho_array: np.ndarray,
    z_mm: float,
    z_sup_mm: float,
    props_sup: OpticalProperties,
    props_cortex: OpticalProperties,
    n_s: int = 500,
    s_max: float = 30.0
) -> np.ndarray:
    """
    Vectorized version of two-layer Kienle fluence for multiple positions at once.

    Computes fluence at (rho_i, z) for all rho_i in rho_array, sharing the same z.
    This is much more efficient when evaluating fluence at many radial positions.

    Parameters:
        rho_array: Array of radial distances (shape: (n_x,))
        z_mm: Depth (scalar)
        z_sup_mm: Thickness of superficial layer
        props_sup: Optical properties of superficial layer
        props_cortex: Optical properties of cortex layer
        n_s: Number of Hankel quadrature points
        s_max: Maximum spatial frequency (mm^-1)

    Returns:
        Fluence values at each rho (shape: (n_x,))
    """
    # Quadrature setup
    s_vals = np.linspace(0.001, s_max, n_s)
    ds = (s_max - 0.001) / (n_s - 1) if n_s > 1 else 1.0

    # Extract optical properties
    D1 = props_sup.D
    D2 = props_cortex.D
    mua1 = props_sup.mua
    mua2 = props_cortex.mua
    musp1 = props_sup.musp

    # Extrapolated boundary
    A = 2.95
    zb = 2.0 * A * D1
    zp = 1.0 / musp1

    # Compute Hankel domain solution for all s (independent of rho)
    alpha1 = np.sqrt((D1 * s_vals**2 + mua1) / D1)
    alpha2 = np.sqrt((D2 * s_vals**2 + mua2) / D2)
    Da = (D1*alpha1 - D2*alpha2) / (D1*alpha1 + D2*alpha2)

    if z_mm < zp:
        arg = zp - z_mm
    else:
        arg = z_mm - zp

    dum1 = np.exp(-alpha1*arg) - np.exp(-alpha1*(2*zb + zp + z_mm))
    dum2 = (np.exp(-alpha1*(-zp + 2*z_sup_mm - z_mm))
            - np.exp(-alpha1*(2*zb + zp + 2*z_sup_mm - z_mm))
            - np.exp(-alpha1*(-zp + 2*z_sup_mm + 2*zb + z_mm))
            + np.exp(-alpha1*(4*zb + zp + 2*z_sup_mm + z_mm)))

    if z_mm < z_sup_mm:
        phi_k = (dum1 + Da*dum2) / (2*D1*alpha1)
    else:
        dum3 = np.exp(alpha2*(z_sup_mm - z_mm)) / (D1*alpha1 + D2*alpha2)
        dum4 = np.exp(alpha1*(zp - z_sup_mm)) - np.exp(-alpha1*(2*zb + zp + z_sup_mm))
        dum5 = (np.exp(alpha1*(zp - 3*z_sup_mm - 2*zb))
                - np.exp(-alpha1*(4*zb + zp + 3*z_sup_mm)))
        phi_k = dum3 * (dum4 - Da*dum5)

    # Vectorized inverse Hankel transform for all rho values
    # j0_vals shape: (n_x, n_s)
    j0_vals = j0(s_vals[None, :] * rho_array[:, None])

    # integrand shape: (n_x, n_s)
    integrand = j0_vals * (phi_k[None, :] * s_vals[None, :])

    # Integrate over s for each rho
    fluence_array = (1.0 / (2.0 * np.pi)) * _np_trapz(integrand, dx=ds, axis=1)

    return np.maximum(fluence_array, 1e-20)


# =============================================================================
# SENSITIVITY / JACOBIAN COMPUTATION
# =============================================================================

def compute_sensitivity_map_two_layer(
    sds_mm: float,
    wavelength_nm: int,
    model: TwoLayerModel,
    z_max_mm: float = 35.0,
    dz_mm: float = 0.5,
    x_range_mm: float = 50.0,
    dx_mm: float = 1.0,
    y_range_mm: float = 50.0,
    dy_mm: float = 2.0,
    n_s: int = 500
) -> Tuple[float, float, float]:
    """
    Compute sensitivity distribution for a source-detector pair using the
    proper two-layer Kienle SDA fluence.

    J(r) ∝ Φ_source(r) × Φ_detector(r)  (adjoint / Born approximation)

    Full 3D volume integral over (x, y, z).
    The integral exploits the symmetry J(x, y, z) = J(x, -y, z) to integrate
    y from 0 to y_max and apply a weight of 2 for all y > 0, halving the
    computational cost relative to the full symmetric range.

    Parameters:
        sds_mm: Source-detector separation (mm)
        wavelength_nm: Wavelength for optical properties
        model: Two-layer tissue model
        z_max_mm: Maximum depth to integrate (mm)
        dz_mm: Depth resolution (mm)
        x_range_mm: Lateral extent to integrate along source-detector axis (mm)
        dx_mm: Lateral resolution along x (mm)
        y_range_mm: Lateral extent to integrate perpendicular to S-D axis (mm)
        dy_mm: Lateral resolution along y (mm); default 2.0 mm for performance
        n_s: Number of Hankel quadrature points

    Returns:
        S_total: Total integrated sensitivity
        S_cortex: Sensitivity in cortex layer only
        f_cortex: Cortical fraction (S_cortex / S_total)
    """
    z_sup = model.z_superficial_mm
    props_sup = model.props_superficial.get(wavelength_nm)
    props_ctx = model.props_cortex.get(wavelength_nm)

    if props_sup is None or props_ctx is None:
        available = list(model.props_superficial.keys())
        nearest = min(available, key=lambda x: abs(x - wavelength_nm))
        props_sup = model.props_superficial[nearest]
        props_ctx = model.props_cortex[nearest]

    print(f"    Computing 3D sensitivity for wavelength {wavelength_nm} nm, "
          f"z_sup={z_sup:.1f} mm...")

    z_vals = np.arange(dz_mm / 2, z_max_mm, dz_mm)
    x_vals = np.arange(-x_range_mm, x_range_mm + dx_mm, dx_mm)
    # This exactly represents the integral over -y_max to +y_max.
    y_vals = np.arange(0.0, y_range_mm + dy_mm / 2, dy_mm)

    S_total = 0.0
    S_cortex = 0.0

    for i_z, z in enumerate(z_vals):
        if i_z % max(1, len(z_vals) // 10) == 0:
            print(f"      z = {z:.1f} mm ({i_z}/{len(z_vals)})")

        for y in y_vals:
            # y=0 lies on the symmetry plane → weight 1; y>0 → weight 2
            y_weight = 1.0 if y == 0.0 else 2.0

            # Radial distances from source (at x=0, y=0) and detector (at x=sds, y=0)
            rho_src_array = np.sqrt(x_vals**2 + y**2)
            rho_det_array = np.sqrt((x_vals - sds_mm)**2 + y**2)

            phi_s_array = fluence_two_layer_kienle_vectorized(
                rho_src_array, z, z_sup, props_sup, props_ctx, n_s=n_s
            )
            phi_d_array = fluence_two_layer_kienle_vectorized(
                rho_det_array, z, z_sup, props_sup, props_ctx, n_s=n_s
            )

            J_array = phi_s_array * phi_d_array
            dV = dx_mm * dy_mm * dz_mm * y_weight

            S_total += np.sum(J_array) * dV
            if z >= z_sup:
                S_cortex += np.sum(J_array) * dV

    f_cortex = S_cortex / S_total if S_total > 0 else 0.0
    print(f"    f_cortex({wavelength_nm} nm) = {f_cortex:.4f}")
    return S_total, S_cortex, float(f_cortex)


def compute_sensitivity_map(
    sds_mm: float,
    wavelength_nm: int,
    model: TwoLayerModel,
    z_max_mm: float = 35.0,
    dz_mm: float = 0.5,
    x_range_mm: float = 50.0,
    dx_mm: float = 1.0,
    y_range_mm: float = 50.0,
    dy_mm: float = 2.0,
    n_s: int = 500
) -> Tuple[float, float, float]:
    """
    Wrapper for sensitivity computation using the two-layer Kienle method.

    Now performs a proper 3D volume integral (see
    compute_sensitivity_map_two_layer for details).
    """
    return compute_sensitivity_map_two_layer(
        sds_mm, wavelength_nm, model,
        z_max_mm, dz_mm,
        x_range_mm, dx_mm,
        y_range_mm, dy_mm,
        n_s
    )


# =============================================================================
# Cortical sensitivity fraction -- SINGLE SOURCE OF TRUTH
# =============================================================================
# f_cortex is read from the versioned production Monte-Carlo file via
# fcortex_source (which loads fcortex_production.json). There is no hard-coded
# table here: the same converged, wavelength-specific values drive both the
# synthetic data generation and the correction, so the two are guaranteed
# consistent. The production run establishes L_max / z_max / annulus / photon
# convergence and reports batch-based uncertainty (see mc_production.py and
# Supplementary Table S1).
import fcortex_source as _fcs


def f_cortex_mc(sds_mm: float, wavelength_nm: int) -> float:
    """Two-layer cortical sensitivity fraction f_cortex(SDS, lambda) from the
    single production Monte-Carlo source (fcortex_production.json)."""
    return _fcs.f_cortex_2L(sds_mm, wavelength_nm)


# =============================================================================
# MBLL INVERSION
# =============================================================================

def mbll_inversion_batch(
    delta_OD_matrix: np.ndarray,
    wavelengths: List[int],
    sds_mm: float,
    DPF: Dict[int, float],
    extinction: ExtinctionCoefficients,
    f_cortex_dict: Optional[Dict[int, float]] = None
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Vectorized MBLL inversion over all time samples at once.

    Parameters:
        delta_OD_matrix: Shape (n_samples, n_wavelengths) OD time series
        wavelengths: List of wavelengths (nm)
        sds_mm: Source-detector separation (mm)
        DPF: Per-wavelength DPF dictionary
        extinction: Extinction coefficient data
        f_cortex_dict: Optional per-wavelength f_cortex (PV correction)

    Returns:
        HbO: (n_samples,) array of ΔHbO₂ estimates
        HbR: (n_samples,) array of ΔHbR estimates
    """
    # Convention note: this implementation works entirely in the DECADIC
    # (base-10) convention -- base-10 optical density with the decadic
    # extinction coefficients -- so there is no explicit ln(10) factor here.
    # The forward model (generate_synthetic_fnirs_data) uses the same decadic
    # convention, so the recovered concentrations are identical to the
    # Napierian form OD = ln(10)*eps*C*L written in the manuscript (Eq. 1):
    # ln(10) simply cancels between the forward model and this inversion.
    E = extinction.get_matrix(wavelengths) * 1e-4      # (n_wl, 2) in mm⁻¹ µM⁻¹
    E_inv = np.linalg.pinv(E)                          # (2, n_wl)
    L_eff = np.array([DPF[wl] * sds_mm for wl in wavelengths])  # (n_wl,)

    OD = delta_OD_matrix.copy()                        # (n_samples, n_wl)

    # Optional PV correction at OD level
    if f_cortex_dict is not None:
        for i, wl in enumerate(wavelengths):
            f = f_cortex_dict.get(wl, 1.0)
            if f > 0.01:
                OD[:, i] /= f

    # Normalise by effective pathlength and invert
    OD_norm = OD / L_eff[None, :]                      # (n_samples, n_wl)
    conc = (E_inv @ OD_norm.T).T                       # (n_samples, 2)

    return conc[:, 0], conc[:, 1]


# =============================================================================
# κ CORRECTION FACTORS
# =============================================================================

def compute_kappa_DPF(
    wavelengths: List[int],
    DPF_used: Dict[int, float],
    DPF_star: Dict[int, float]
) -> float:
    """Compute κ(DPF) - pathlength calibration correction."""
    dpf_used_mean = np.mean([DPF_used.get(w, 6.0) for w in wavelengths])
    dpf_star_mean = np.mean([DPF_star.get(w, 6.0) for w in wavelengths])
    if dpf_star_mean <= 0:
        return 1.0
    return dpf_used_mean / dpf_star_mean


def compute_kappa_PV(
    sds_mm: float,
    wavelengths: List[int],
    model: TwoLayerModel,
    n_s: int = 500
) -> Tuple[Dict[int, float], float]:
    """
    Compute wavelength-specific κ(PV) - partial volume correction.

    Returns per-wavelength f_cortex values and mean κ(PV).

    Returns:
        f_cortex_dict: Dict mapping wavelength → f_cortex
        kappa_PV_mean: Mean correction factor across wavelengths
    """
    f_cortex_dict = {}
    for wl in wavelengths:
        f_ctx = f_cortex_mc(sds_mm, wl)  # stable MC f_cortex (was non-converged Hankel)
        f_cortex_dict[wl] = f_ctx

    f_cortex_mean = np.mean(list(f_cortex_dict.values()))
    if f_cortex_mean <= 0.01:
        kappa_PV_mean = 100.0
    else:
        kappa_PV_mean = 1.0 / f_cortex_mean

    return f_cortex_dict, kappa_PV_mean


def _compute_R2_SS_single(
    y_long: np.ndarray,
    y_short: np.ndarray,
    task_regressor: Optional[np.ndarray],
    sfreq: float,
    band: Tuple[float, float] = (0.01, 0.5)
) -> float:
    """
    Ordinary coefficient of determination (R^2) from the SSR regression of the
    long-channel dOD on the short-channel dOD, i.e. the fraction of long-channel
    variance removed by SSR:

        R2_SS = 1 - Var(residual) / Var(dOD_long),

    where residual = dOD_long - [1, dOD_short] @ beta (OLS).  This is the SAME
    statistic the real-data pipeline uses (fnirs_kappa_realdata_v2.py) and the
    definition stated in the manuscript methods, so the synthetic and real-data
    V_SSR diagnostics are now harmonized.

    DEPRECATED / NOT USED FOR THE REPORTED DIAGNOSTIC.  This helper band-pass
    filters the signals to `band` (default 0.01-0.5 Hz) before regressing, so its
    R^2 is a band-limited long-short coupling statistic, NOT the R^2 of the actual
    SSR operator.  The reported synthetic V_SSR is now built directly from the R^2
    returned by perform_ssr() (unfiltered OLS long-on-short), matching the
    real-data pipeline and the manuscript definition.  This function is retained
    only for reference/backward compatibility.

    Note: this is deliberately NOT a task-partialled (partial) R^2.  The task
    regressor is accepted for call-signature compatibility but is not used.
    Works on any 1-D signal pair (OD at a single wavelength, or concentration).
    Returns R2_SS clipped to [0, 0.99].
    """
    del task_regressor  # harmonized ordinary R^2 does not partial out the task

    nyq = sfreq / 2
    if band[1] >= nyq:
        band = (band[0], nyq * 0.95)

    y_long_f = bandpass_filter(y_long, sfreq, band[0], band[1])

    if y_short.ndim == 1:
        y_short = y_short.reshape(-1, 1)

    y_short_f = np.apply_along_axis(
        lambda x: bandpass_filter(x, sfreq, band[0], band[1]),
        0, y_short
    )

    n = len(y_long_f)
    # Ordinary OLS regression of the long channel on the short channel(s) with
    # an intercept -- exactly the SSR operator, not a task-conditioned model.
    X = np.column_stack([np.ones(n), y_short_f])
    beta = pinv(X) @ y_long_f
    resid = y_long_f - X @ beta

    SS_res = float(np.sum(resid**2))
    SS_tot = float(np.sum((y_long_f - np.mean(y_long_f))**2))

    if SS_tot > 0:
        R2_SS = 1.0 - SS_res / SS_tot
        R2_SS = float(np.clip(R2_SS, 0, 0.99))
    else:
        R2_SS = 0.0

    return float(R2_SS)


def compute_V_SSR(
    long_OD_dict: Dict[int, np.ndarray],
    short_OD_dict: Dict[int, np.ndarray],
    wavelengths: List[int],
    task_regressor: Optional[np.ndarray],
    sfreq: float,
    band: Tuple[float, float] = (0.01, 0.5)
) -> Tuple[Dict[int, float], Dict[int, float], float, float]:
    """
    DEPRECATED -- band-limited diagnostic only; DO NOT USE for the reported V_SSR.

    This function calls the band-pass-filtered helper _compute_R2_SS_single()
    (default 0.01-0.5 Hz), so its R² is a band-limited long-short coupling
    statistic, NOT the R² of the actual SSR operator.  The reported synthetic
    V_SSR is built directly from the R² that perform_ssr() returns (unfiltered
    OLS long-on-short with an intercept), which is the exact fraction of
    long-channel variance removed by the SSR operation the pipeline applies.
    This function is retained only for reference/backward compatibility; prefer
    perform_ssr()'s returned R² for any reported diagnostic.

    Compute V(SSR) per wavelength at the OD level.

    V_SSR(λ) = 1 / (1 - R²_SS(λ))

    where R²_SS(λ) is the (band-limited) ordinary coefficient of determination
    from the OLS regression of the *uncorrected* long-channel ΔOD(λ) on the
    short-channel ΔOD(λ) at the same wavelength.
    It is NOT a task-partialled R², and it must NOT be computed against the
    SSR-corrected long-channel signal — by construction of OLS, the residual of
    the SSR fit is orthogonal to the short-channel regressor, so doing so would
    force R²_SS ≈ 0 and V_SSR ≈ 1 regardless of the true coupled variance.

    Because the absorption and scattering coefficients differ at each
    wavelength, the superficial contamination fraction — and hence the
    SSR regression R² — is wavelength-dependent.  Computing V_SSR(λ)
    per wavelength preserves this physical dependence and is consistent
    with the per-wavelength OD-level SSR and PV corrections.

    V_SSR is reported PER WAVELENGTH only.  A cross-wavelength arithmetic
    mean is returned (last two outputs) for backward compatibility, but it is
    NOT used in the results tables, figures or summaries, because the
    per-wavelength values differ substantially (e.g. 760 vs 850 nm) so their
    mean is not a physically meaningful quantity.

    Parameters:
        long_OD_dict: *Uncorrected* long-channel ΔOD(λ) per wavelength
        short_OD_dict: Short-channel ΔOD(λ) per wavelength
        wavelengths: List of wavelengths (nm)
        task_regressor: Convolved task design regressor (or None)
        sfreq: Sampling frequency (Hz)
        band: Bandpass filter range (Hz)

    Returns:
        V_SSR_dict: {wavelength: V_SSR(λ)}
        R2_SS_dict: {wavelength: R²_SS(λ)}
        V_SSR_mean: Arithmetic mean of V_SSR across wavelengths
        R2_SS_mean: Arithmetic mean of R²_SS across wavelengths
    """
    V_SSR_dict: Dict[int, float] = {}
    R2_SS_dict: Dict[int, float] = {}

    for wl in wavelengths:
        R2 = _compute_R2_SS_single(
            long_OD_dict[wl], short_OD_dict[wl],
            task_regressor, sfreq, band
        )
        R2_SS_dict[wl] = R2
        V_SSR_dict[wl] = 1.0 / (1.0 - R2)

    V_SSR_mean = float(np.mean(list(V_SSR_dict.values())))
    R2_SS_mean = float(np.mean(list(R2_SS_dict.values())))

    return V_SSR_dict, R2_SS_dict, V_SSR_mean, R2_SS_mean


# =============================================================================
# SYNTHETIC DATA GENERATION
# =============================================================================

def generate_hrf(t: np.ndarray, peak_time: float = 6.0, undershoot: float = 0.35) -> np.ndarray:
    """Generate a canonical hemodynamic response function (HRF)."""
    a1, b1 = peak_time, 1.0
    a2, b2 = 16.0, 1.0

    hrf = (t**(a1-1) * np.exp(-t/b1) / (b1**a1 * math.gamma(a1)) -
           undershoot * t**(a2-1) * np.exp(-t/b2) / (b2**a2 * math.gamma(a2)))
    hrf = hrf / np.max(np.abs(hrf))
    return hrf


def generate_synthetic_fnirs_data(
    n_subjects: int = 5,
    duration_s: float = 660.0,
    sfreq: float = 7.8125,
    sds_long_mm: List[float] = [25.0, 30.0, 35.0, 38.0, 40.0],
    sds_short_mm: float = 10.0,
    wavelengths: List[int] = [760, 850],
    model: Optional[TwoLayerModel] = None,
    seed: int = 42
) -> Dict:
    """
    Generate synthetic fNIRS data with known ground truth.

    The synthetic paradigm is LOOSELY INSPIRED BY the BIDS-NIRS-Tapping dataset
    (Luke et al., 2021) used in the real-data proof-of-concept, but it intentionally
    uses 60 regular 5-s-on / 5-s-off blocks for controlled validation. It is NOT a
    temporal reproduction of the experimental paradigm: the real dataset has three
    conditions, ~90 trials, and irregular inter-onset intervals, which this idealized
    regular block design does not reproduce. Shared parameters are 5 subjects,
    sampling rate 7.8125 Hz, and wavelengths 760 and 850 nm. The design is used for
    qualitative pipeline comparison rather than direct quantitative design matching.

    The synthetic data includes:
    1. True cortical activation (HbO ↑, HbR ↓)
    2. Systemic/superficial fluctuations (Mayer waves, respiration, cardiac)
    3. Measurement noise

    Each wavelength's ΔOD is generated using the wavelength-specific
    cortical fraction f_cortex(λ), so the forward model correctly exhibits
    the wavelength-dependent partial-volume effect that the correction removes.
    Noise is added at the OD level as a fixed, separation-independent additive
    Gaussian term (assumed additive OD measurement noise, NOT a photon-count /
    shot-noise model, whose variance would grow with source-detector separation).

    Short-channel ΔOD is stored per wavelength so that
    the SSR step can regress each long-channel ΔOD(λ) against the
    short-channel ΔOD at the same wavelength.

    Returns:
        Dictionary containing all synthetic data and ground truth values
    """
    np.random.seed(seed)

    if model is None:
        model = TwoLayerModel.default_adult(z_superficial_mm=12.0)

    extinction = ExtinctionCoefficients()
    n_samples = int(duration_s * sfreq)
    t = np.arange(n_samples) / sfreq

    # Task design: block design with 5s on, 5s off
    # (idealized regular block design loosely inspired by BIDS-NIRS-Tapping: 60 blocks
    #  of 5 s task / 5 s rest; NOT a temporal reproduction of the real irregular paradigm)
    block_duration = 5.0
    task_onsets = np.arange(30, duration_s - 30, 2 * block_duration)

    # Generate HRF
    hrf_t = np.arange(0, 30, 1/sfreq)
    hrf = generate_hrf(hrf_t)

    # Task regressor (convolved with HRF)
    task_boxcar = np.zeros(n_samples)
    for onset in task_onsets:
        start = int(onset * sfreq)
        end = int((onset + block_duration) * sfreq)
        task_boxcar[start:min(end, n_samples)] = 1.0

    task_regressor = np.convolve(task_boxcar, hrf, mode='full')[:n_samples]
    task_regressor = task_regressor / np.max(task_regressor)

    # Storage
    data = {
        'n_subjects': n_subjects,
        'n_samples': n_samples,
        'sfreq': sfreq,
        't': t,
        'task_regressor': task_regressor,
        'wavelengths': wavelengths,
        'sds_long_mm': sds_long_mm,
        'sds_short_mm': sds_short_mm,
        'model': model,
        'subjects': []
    }

    # True cortical activation amplitude (µM)
    true_HbO_cortex = 2.5
    true_HbR_cortex = -0.8

    dpf_true = get_dpf_dict(wavelengths)

    print("  Generating synthetic data...")
    for subj in range(n_subjects):
        subj_data = {
            'id': f'sub-{subj+1:02d}',
            'channels': []
        }

        # Subject-specific variations
        amplitude_scale = 0.8 + 0.4 * np.random.rand()

        # Systemic fluctuations (same for all channels within a subject)
        # Use amplitude-modulated physiological oscillations to produce
        # realistic channel correlations after bandpass filtering.
        am_mayer   = 1 + 0.4 * np.sin(2*np.pi*0.005*t + 2*np.pi*np.random.rand())
        am_resp    = 1 + 0.3 * np.sin(2*np.pi*0.007*t + 2*np.pi*np.random.rand())
        am_cardiac = 1 + 0.3 * np.sin(2*np.pi*0.013*t + 2*np.pi*np.random.rand())

        mayer   = 0.3 * am_mayer   * np.sin(2*np.pi*0.10*t + 2*np.pi*np.random.rand())
        resp    = 0.2 * am_resp    * np.sin(2*np.pi*0.25*t + 2*np.pi*np.random.rand())
        cardiac = 0.1 * am_cardiac * np.sin(2*np.pi*1.00*t + 2*np.pi*np.random.rand())
        # Additional VLF autonomic fluctuation (0.02-0.05 Hz)
        vlf     = 0.2 * np.sin(2*np.pi*0.04*t + 2*np.pi*np.random.rand())

        systemic_HbO = mayer + resp + cardiac + vlf
        systemic_HbR = -0.3 * (mayer + resp + cardiac + vlf)

        # ----------------------------------------------------------------
        # Short channel (superficial signal only)
        # ----------------------------------------------------------------
        short_HbO = systemic_HbO.copy()
        short_HbR = systemic_HbR.copy()

        # store per-wavelength ΔOD for the short channel
        # so the SSR regression can operate on the same quantity (OD) as the
        # long channel, at the same wavelength.
        # Noise is added at the OD level for consistency with the long channels.
        noise_level_OD_short = 0.001  # same σ_OD as long channels
        short_delta_OD: Dict[int, np.ndarray] = {}
        for wl in wavelengths:
            eps_wl = extinction.coefficients[wl]
            eps_HbO2_mm_wl = eps_wl[0] * 1e-4
            eps_HbR_mm_wl  = eps_wl[1] * 1e-4
            L_eff_short = dpf_true[wl] * sds_short_mm
            clean_OD_short = (
                eps_HbO2_mm_wl * short_HbO + eps_HbR_mm_wl * short_HbR
            ) * L_eff_short
            short_delta_OD[wl] = clean_OD_short + noise_level_OD_short * np.random.randn(n_samples)

        subj_data['short_channel'] = {
            'sds_mm': sds_short_mm,
            'HbO': short_HbO,
            'HbR': short_HbR,
            'delta_OD': short_delta_OD,
            'truth_HbO_superficial': systemic_HbO,
            'truth_HbR_superficial': systemic_HbR
        }

        # ----------------------------------------------------------------
        # Long channels
        # ----------------------------------------------------------------
        for sds in sds_long_mm:
            # Compute wavelength-specific f_cortex for this SDS
            f_cortex_dict: Dict[int, float] = {}
            for wl in wavelengths:
                f_ctx = f_cortex_mc(sds, wl)  # stable MC f_cortex (was non-converged Hankel)
                f_cortex_dict[wl] = f_ctx

            # Mean f_cortex (used only for reporting / stored Hb traces)
            f_cortex_mean = float(np.mean(list(f_cortex_dict.values())))

            # True cortical time series (same for all wavelengths)
            cortex_HbO = amplitude_scale * true_HbO_cortex * task_regressor
            cortex_HbR = amplitude_scale * true_HbR_cortex * task_regressor
            # f_cortex(λ).  This ensures the synthetic data contains the
            # actual wavelength-dependent partial-volume bias that the
            # correction method is designed to remove.
            delta_OD: Dict[int, np.ndarray] = {}
            noise_level_OD = 0.001  # σ = 0.001 OD, assumed additive OD measurement noise (fixed, separation-independent)

            for wl in wavelengths:
                f_ctx_wl  = f_cortex_dict[wl]
                f_sup_wl  = 1.0 - f_ctx_wl

                # Wavelength-specific mixed Hb concentrations (µM)
                HbO_wl = f_ctx_wl * cortex_HbO + f_sup_wl * systemic_HbO
                HbR_wl = f_ctx_wl * cortex_HbR + f_sup_wl * systemic_HbR

                eps_wl = extinction.coefficients[wl]
                eps_HbO2_mm = eps_wl[0] * 1e-4
                eps_HbR_mm  = eps_wl[1] * 1e-4
                L_eff = dpf_true[wl] * sds

                # Compute clean OD, then add noise at the OD level
                # (assumed additive OD measurement noise: fixed, separation-independent)
                clean_OD = (eps_HbO2_mm * HbO_wl + eps_HbR_mm * HbR_wl) * L_eff
                delta_OD[wl] = clean_OD + noise_level_OD * np.random.randn(n_samples)

            # Measured Hb traces for storage/display (mean f_cortex approx.)
            # Noise is added at OD level above; these traces approximate the
            # noisy signal by back-converting from the mean-wavelength OD noise
            f_sup_mean = 1.0 - f_cortex_mean
            measured_HbO = (f_cortex_mean * cortex_HbO
                            + f_sup_mean * systemic_HbO)
            measured_HbR = (f_cortex_mean * cortex_HbR
                            + f_sup_mean * systemic_HbR)

            channel_data = {
                'sds_mm': sds,
                'f_cortex': f_cortex_mean,
                'f_cortex_dict': f_cortex_dict,
                'delta_OD': delta_OD,
                'truth_HbO_cortex': cortex_HbO,
                'truth_HbR_cortex': cortex_HbR,
                'truth_HbO_superficial': systemic_HbO,
                'truth_HbR_superficial': systemic_HbR,
                'measured_HbO': measured_HbO,
                'measured_HbR': measured_HbR
            }
            subj_data['channels'].append(channel_data)

        data['subjects'].append(subj_data)

    data['true_HbO_cortex'] = true_HbO_cortex
    data['true_HbR_cortex'] = true_HbR_cortex

    return data


# =============================================================================
# SHORT-SEPARATION REGRESSION
# =============================================================================

def perform_ssr(y_long: np.ndarray, y_short: np.ndarray) -> Tuple[np.ndarray, float]:
    """
    Perform Short-Separation Regression (SSR) to remove systemic contamination.

    Parameters:
        y_long: Long-separation channel time series
        y_short: Short-separation channel time series

    Returns:
        y_corrected: SSR-corrected long channel signal
        R2: Variance explained by short channel
    """
    n = len(y_long)
    X = np.column_stack([np.ones(n), y_short])
    beta = pinv(X) @ y_long
    y_pred = X @ beta
    y_corrected = y_long - beta[1] * y_short

    SS_total = np.sum((y_long - np.mean(y_long))**2)
    SS_resid = np.sum((y_long - y_pred)**2)
    R2 = 1 - SS_resid / SS_total if SS_total > 0 else 0

    return y_corrected, R2


# =============================================================================
# ANALYSIS PIPELINE
# =============================================================================

def run_analysis_pipeline(data: Dict, n_s: int = 500) -> pd.DataFrame:
    """
    Run the complete analysis pipeline on synthetic data with wavelength-specific
    partial volume correction at the OD level.

    Pipeline:
    1. MBLL inversion (biased, no correction; kept for comparison)
    2. SSR at the OD level for each wavelength using the short-channel ΔOD at
       the same wavelength λ
    3. Apply wavelength-specific f_cortex correction: ΔOD_ssr / f_cortex(λ)
    4. MBLL inversion on corrected ODs → corrected concentrations
    5. Compute and report κ factors:
         κ(PV) = 1/f_cortex(λ) is applied at the OD level
         V(SSR, λ) = 1/(1−R²_SS(λ)) is reported PER WAVELENGTH as a diagnostic
         only (not applied); no cross-wavelength mean or κ_total is formed

    Parameters:
        data: Synthetic data dictionary
        n_s: Number of Hankel quadrature points

    Returns:
        DataFrame with results for all subjects and channels
    """
    print("\n  Running analysis pipeline...")
    extinction = ExtinctionCoefficients()
    wavelengths = data['wavelengths']
    model = data['model']
    sfreq = data['sfreq']
    task_regressor = data['task_regressor']

    results = []
    DPF_wavelength = get_dpf_dict(wavelengths)
    DPF_used = DPF_wavelength
    DPF_star = DPF_wavelength

    for subj_data in data['subjects']:
        subj_id = subj_data['id']
        short_ch = subj_data['short_channel']

        for ch_data in subj_data['channels']:
            sds = ch_data['sds_mm']
            f_cortex_true = ch_data['f_cortex']
            f_cortex_dict = ch_data['f_cortex_dict']

            truth_HbO = ch_data['truth_HbO_cortex']
            truth_HbR = ch_data['truth_HbR_cortex']

            X_task = np.column_stack([np.ones(len(task_regressor)), task_regressor])
            beta_truth_HbO = pinv(X_task) @ truth_HbO
            beta_truth_HbR = pinv(X_task) @ truth_HbR
            truth_beta_HbO = beta_truth_HbO[1]
            truth_beta_HbR = beta_truth_HbR[1]

            # ------------------------------------------------------------------
            # Step 1: Standard MBLL inversion (no PV correction) — vectorized
            # ------------------------------------------------------------------
            delta_OD_matrix = np.column_stack(
                [ch_data['delta_OD'][wl] for wl in wavelengths]
            )
            mbll_HbO, mbll_HbR = mbll_inversion_batch(
                delta_OD_matrix, wavelengths, sds, DPF_wavelength,
                extinction, f_cortex_dict=None
            )

            beta_mbll_HbO = pinv(X_task) @ mbll_HbO
            beta_mbll_HbR = pinv(X_task) @ mbll_HbR
            mbll_beta_HbO = beta_mbll_HbO[1]
            mbll_beta_HbR = beta_mbll_HbR[1]

            # ------------------------------------------------------------------
            # Step 2: SSR at OD level
            # at the SAME wavelength, preserving Beer-Lambert linearity.
            # ------------------------------------------------------------------
            ssr_OD_dict: Dict[int, np.ndarray] = {}
            ssr_R2_dict: Dict[int, float] = {}
            for wl in wavelengths:
                short_OD_wl = short_ch['delta_OD'][wl]   # same-wavelength OD
                ssr_OD, ssr_r2 = perform_ssr(ch_data['delta_OD'][wl], short_OD_wl)
                ssr_OD_dict[wl] = ssr_OD
                # Capture the R^2 of the ACTUAL SSR regression (unfiltered OD,
                # OLS with intercept) so the reported V_SSR is exactly the
                # fraction of long-channel variance removed by this operation.
                ssr_R2_dict[wl] = float(np.clip(ssr_r2, 0.0, 0.99))

            # ------------------------------------------------------------------
            # Step 3: Wavelength-specific PV correction at OD level
            #         ΔOD_corrected(λ) = ΔOD_ssr(λ) / f_cortex(λ)
            # ------------------------------------------------------------------
            corrected_OD_dict: Dict[int, np.ndarray] = {}
            for wl in wavelengths:
                f_cortex_wl = f_cortex_dict.get(wl, 1.0)
                if f_cortex_wl > 0.01:
                    corrected_OD_dict[wl] = ssr_OD_dict[wl] / f_cortex_wl
                else:
                    corrected_OD_dict[wl] = ssr_OD_dict[wl]

            # ------------------------------------------------------------------
            # Step 4: MBLL inversion on corrected OD — vectorized
            # ------------------------------------------------------------------
            corrected_OD_matrix = np.column_stack(
                [corrected_OD_dict[wl] for wl in wavelengths]
            )
            corrected_HbO, corrected_HbR = mbll_inversion_batch(
                corrected_OD_matrix, wavelengths, sds, DPF_wavelength,
                extinction, f_cortex_dict=None
            )

            beta_corrected_HbO = pinv(X_task) @ corrected_HbO
            beta_corrected_HbR = pinv(X_task) @ corrected_HbR
            corrected_beta_HbO = beta_corrected_HbO[1]
            corrected_beta_HbR = beta_corrected_HbR[1]

            # ------------------------------------------------------------------
            # Step 5: Compute κ factors
            # ------------------------------------------------------------------
            kappa_DPF = compute_kappa_DPF(wavelengths, DPF_used, DPF_star)

            f_cortex_computed_dict, kappa_PV = compute_kappa_PV(
                sds, wavelengths, model, n_s=n_s
            )
            f_cortex_computed = float(
                np.mean(list(f_cortex_computed_dict.values()))
            )
            # V_SSR is built DIRECTLY from the R^2 returned by the actual SSR
            # operation (perform_ssr above): R2_SS(λ) is the ordinary
            # coefficient of determination of the same unfiltered OLS
            # long-on-short regression (with intercept) that produced the
            # SSR-corrected signal, so V_SSR(λ) = 1/(1-R2_SS(λ)) is exactly the
            # fraction of long-channel variance removed by SSR at that
            # wavelength.  This is the SAME statistic the real-data pipeline
            # reports.  It is a per-wavelength DIAGNOSTIC and is NOT applied
            # (only kappa_PV is applied at the OD level); no cross-wavelength
            # mean or composite kappa_total is formed because the 760/850 nm
            # values differ substantially and their mean is not meaningful.
            R2_SS_dict = {wl: ssr_R2_dict[wl] for wl in wavelengths}
            V_SSR_dict = {wl: 1.0 / (1.0 - R2_SS_dict[wl]) for wl in wavelengths}

            # ------------------------------------------------------------------
            # Compute errors
            # ------------------------------------------------------------------
            error_mbll_HbO = mbll_beta_HbO - truth_beta_HbO
            error_mbll_HbR = mbll_beta_HbR - truth_beta_HbR
            error_corr_HbO = corrected_beta_HbO - truth_beta_HbO
            error_corr_HbR = corrected_beta_HbR - truth_beta_HbR

            corr_short, _ = pearsonr(mbll_HbO, short_ch['HbO'])

            # Build per-wavelength V_SSR columns
            V_SSR_per_wl = {
                f'V_SSR_{wl}': V_SSR_dict[wl] for wl in wavelengths
            }
            R2_SS_per_wl = {
                f'R2_SS_{wl}': R2_SS_dict[wl] for wl in wavelengths
            }

            results.append({
                'subject': subj_id,
                'sds_mm': sds,
                'f_cortex_true': f_cortex_true,
                'f_cortex_computed': f_cortex_computed,
                'truth_beta_HbO': truth_beta_HbO,
                'truth_beta_HbR': truth_beta_HbR,
                'mbll_beta_HbO': mbll_beta_HbO,
                'mbll_beta_HbR': mbll_beta_HbR,
                'corrected_beta_HbO': corrected_beta_HbO,
                'corrected_beta_HbR': corrected_beta_HbR,
                'kappa_DPF': kappa_DPF,
                'kappa_PV': kappa_PV,
                **V_SSR_per_wl,
                **R2_SS_per_wl,
                'error_mbll_HbO': error_mbll_HbO,
                'error_mbll_HbR': error_mbll_HbR,
                'error_corr_HbO': error_corr_HbO,
                'error_corr_HbR': error_corr_HbR,
                'corr_short': corr_short
            })

    return pd.DataFrame(results)


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_results(df: pd.DataFrame, save_path: Optional[str] = None):
    """Create comprehensive 6-panel visualization of results."""
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))

    sds_unique = sorted(df['sds_mm'].unique())

    # 1. Cortical fraction vs SDS
    ax = axes[0, 0]
    f_cortex_mean = [df[df['sds_mm']==s]['f_cortex_computed'].mean() for s in sds_unique]
    ax.plot(sds_unique, f_cortex_mean, 'bo-', markersize=10, linewidth=2)
    ax.set_xlabel('Source-Detector Separation (mm)', fontsize=12)
    ax.set_ylabel('Cortical Fraction (f_cortex)', fontsize=12)
    ax.set_title('(A) Cortical Sensitivity vs SDS', fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.set_ylim([0, 1])

    # 2. κ(PV) vs SDS
    ax = axes[0, 1]
    for sds in sds_unique:
        subset = df[df['sds_mm'] == sds]
        ax.scatter([sds]*len(subset), subset['kappa_PV'], alpha=0.5,
                   label='κ(PV)' if sds == sds_unique[0] else '')
    ax.set_xlabel('Source-Detector Separation (mm)', fontsize=12)
    ax.set_ylabel('κ(PV)', fontsize=12)
    ax.set_title('(B) Partial Volume Correction Factor', fontsize=14)
    ax.grid(True, alpha=0.3)

    # 3. Estimation error vs SDS
    ax = axes[0, 2]
    for sds in sds_unique:
        subset = df[df['sds_mm'] == sds]
        ax.scatter([sds]*len(subset), subset['error_mbll_HbO'], alpha=0.5, c='red',
                   label='MBLL' if sds == sds_unique[0] else '')
        ax.scatter([sds]*len(subset), subset['error_corr_HbO'], alpha=0.5, c='green',
                   label='Corrected' if sds == sds_unique[0] else '')
    ax.axhline(0, color='k', linestyle='--', alpha=0.5)
    ax.set_xlabel('Source-Detector Separation (mm)', fontsize=12)
    ax.set_ylabel('HbO Amplitude Error (µM)', fontsize=12)
    ax.set_title('(C) Estimation Error: MBLL vs Corrected', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 4. True vs Estimated
    ax = axes[1, 0]
    ax.scatter(df['truth_beta_HbO'], df['mbll_beta_HbO'], alpha=0.5, c='red', label='MBLL')
    ax.scatter(df['truth_beta_HbO'], df['corrected_beta_HbO'], alpha=0.5, c='green', label='Corrected')
    lim = [df['truth_beta_HbO'].min() - 0.5, df['truth_beta_HbO'].max() + 0.5]
    ax.plot(lim, lim, 'k--', alpha=0.5, label='Identity')
    ax.set_xlabel('True HbO Amplitude (µM)', fontsize=12)
    ax.set_ylabel('Estimated HbO Amplitude (µM)', fontsize=12)
    ax.set_title('(D) True vs Estimated HbO', fontsize=14)
    ax.legend()
    ax.grid(True, alpha=0.3)

    # 5. κ components.  V_SSR is shown PER WAVELENGTH (diagnostic, not applied);
    #    no cross-wavelength mean is plotted (the 760/850 nm values differ too
    #    much for their average to be meaningful).
    ax = axes[1, 1]
    kappa_cols = ['kappa_DPF', 'kappa_PV', 'V_SSR_760', 'V_SSR_850']
    kappa_means = df.groupby('sds_mm')[kappa_cols].mean()
    x = np.arange(len(sds_unique))
    width = 0.2
    ax.bar(x - 1.5*width, kappa_means['kappa_DPF'], width, label='κ(DPF)', alpha=0.8)
    ax.bar(x - 0.5*width, kappa_means['kappa_PV'], width, label='κ(PV)', alpha=0.8)
    ax.bar(x + 0.5*width, kappa_means['V_SSR_760'], width, label='V$_{SSR}$(760)', alpha=0.8)
    ax.bar(x + 1.5*width, kappa_means['V_SSR_850'], width, label='V$_{SSR}$(850)', alpha=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{s}mm' for s in sds_unique])
    ax.set_xlabel('Source-Detector Separation', fontsize=12)
    ax.set_ylabel('Dimensionless factor / diagnostic', fontsize=12)
    ax.set_title('(E) Applied factors ($\\kappa$) and V$_{SSR}$ diagnostic by SDS', fontsize=14)
    ax.legend()
    ax.axhline(1, color='gray', linestyle='--', alpha=0.5)
    ax.grid(True, alpha=0.3, axis='y')

    # 6. Error reduction summary
    ax = axes[1, 2]
    rmse_mbll = np.sqrt((df['error_mbll_HbO']**2).mean())
    rmse_corr = np.sqrt((df['error_corr_HbO']**2).mean())
    improvement = (1 - rmse_corr/rmse_mbll) * 100

    bars = ax.bar(['MBLL', 'Corrected'], [rmse_mbll, rmse_corr],
                  color=['red', 'green'], alpha=0.7)
    ax.set_ylabel('RMSE (µM)', fontsize=12)
    ax.set_title(f'(F) Error Reduction: {improvement:.1f}%', fontsize=14)
    ax.grid(True, alpha=0.3, axis='y')

    for bar, val in zip(bars, [rmse_mbll, rmse_corr]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.05,
                f'{val:.3f}', ha='center', fontsize=11)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Figure saved to: {save_path}")
    plt.close()
    return fig


def plot_time_series_example(data: Dict, df: pd.DataFrame, subject_idx: int = 0,
                            channel_idx: int = 1, save_path: Optional[str] = None):
    """Plot example time series showing raw MBLL vs corrected."""
    subj = data['subjects'][subject_idx]
    ch = subj['channels'][channel_idx]
    t = data['t']
    task = data['task_regressor']

    sds = ch['sds_mm']
    subj_id = subj['id']
    row = df[(df['subject'] == subj_id) & (df['sds_mm'] == sds)].iloc[0]

    fig, axes = plt.subplots(4, 1, figsize=(14, 13), sharex=True)

    # 1. Task design
    ax = axes[0]
    ax.fill_between(t, 0, task, alpha=0.3, label='Task regressor')
    ax.set_ylabel('Task', fontsize=12)
    ax.set_title(f'{subj_id}, SDS = {sds}mm, f_cortex = {row["f_cortex_computed"]:.2f}', fontsize=14)
    ax.legend(loc='upper right')
    ax.set_ylim([-0.1, 1.2])

    # --- Recompute the MBLL and corrected time series (both HbO2 and HbR) ---
    extinction = ExtinctionCoefficients()
    wavelengths = data['wavelengths']
    DPF_wl = get_dpf_dict(wavelengths)

    # Raw MBLL time series (no correction) -> the "biased" traces.
    delta_OD_matrix = np.column_stack([ch['delta_OD'][wl] for wl in wavelengths])
    mbll_HbO, mbll_HbR = mbll_inversion_batch(delta_OD_matrix, wavelengths, sds, DPF_wl, extinction)

    # Corrected time series, computed with the SAME pipeline as run_analysis_pipeline:
    # per-wavelength SSR at the OD level -> divide by f_cortex(lambda) -> MBLL.
    # We plot these ACTUAL corrected time series (not beta_corrected x task), so every
    # trace is a real time series shown on equal footing.  The corrected traces are
    # therefore as noisy as the MBLL traces, only larger in amplitude (and with the
    # superficial component regressed out) -- they do not impose the task shape.
    short_ch = subj['short_channel']
    f_cortex_dict = {wl: f_cortex_mc(sds, wl) for wl in wavelengths}
    corrected_OD_cols = []
    for wl in wavelengths:
        ssr_OD, _ = perform_ssr(ch['delta_OD'][wl], short_ch['delta_OD'][wl])
        f_ctx = f_cortex_dict[wl]
        corrected_OD_cols.append(ssr_OD / f_ctx if f_ctx > 0.01 else ssr_OD)
    corrected_OD_matrix = np.column_stack(corrected_OD_cols)
    corrected_HbO, corrected_HbR = mbll_inversion_batch(
        corrected_OD_matrix, wavelengths, sds, DPF_wl, extinction, f_cortex_dict=None)

    truth_HbO = ch['truth_HbO_cortex']
    truth_HbR = ch['truth_HbR_cortex']

    # 2. HbO2 comparison
    ax = axes[1]
    ax.plot(t, truth_HbO, 'k-', linewidth=2, label='Truth (cortical)', alpha=0.8)
    ax.plot(t, mbll_HbO, 'r-', linewidth=1.5, label='MBLL (biased)', alpha=0.7)
    ax.plot(t, corrected_HbO, 'g-', linewidth=1.5, label='Corrected (SSR+PV)', alpha=0.7)
    ax.set_ylabel('ΔHbO₂ (µM)', fontsize=12)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 3. HbR comparison (same treatment as HbO2: real time series, not beta x task)
    ax = axes[2]
    ax.plot(t, truth_HbR, 'k-', linewidth=2, label='Truth (cortical)', alpha=0.8)
    ax.plot(t, mbll_HbR, 'r-', linewidth=1.5, label='MBLL (biased)', alpha=0.7)
    ax.plot(t, corrected_HbR, 'g-', linewidth=1.5, label='Corrected (SSR+PV)', alpha=0.7)
    ax.set_ylabel('ΔHbR (µM)', fontsize=12)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    # 4. Short channel
    ax = axes[3]
    ax.plot(t, subj['short_channel']['HbO'], 'b-', linewidth=1,
            label='Short channel (superficial)', alpha=0.7)
    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('ΔHbO₂ (µM)', fontsize=12)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"  Figure saved to: {save_path}")
    plt.close()
    return fig


def plot_robustness(thickness_results: Dict, opt_results: Dict,
                    save_path: Optional[str] = None):
    """Create 4-panel robustness/sensitivity analysis figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # Panel A: f_cortex vs superficial thickness
    ax = axes[0, 0]
    thicknesses_num = sorted([t for t in thickness_results.keys() if isinstance(t, (int, float))])
    f_cortex_vals = [thickness_results[t]['f_cortex'] for t in thicknesses_num]
    ax.plot(thicknesses_num, f_cortex_vals, 'bo-', linewidth=2, markersize=8)
    ax.set_xlabel('Superficial Layer Thickness (mm)', fontsize=12)
    ax.set_ylabel('Cortical Fraction (f_cortex)', fontsize=12)
    ax.set_title('(A) Sensitivity to Superficial Thickness', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Panel B: κ(PV) vs superficial thickness
    ax = axes[0, 1]
    kappa_vals = [thickness_results[t]['kappa_pv'] for t in thicknesses_num]
    ax.plot(thicknesses_num, kappa_vals, 'rs-', linewidth=2, markersize=8)
    ax.set_xlabel('Superficial Layer Thickness (mm)', fontsize=12)
    ax.set_ylabel('Partial Volume Correction (κ(PV))', fontsize=12)
    ax.set_title('(B) κ(PV) vs Superficial Thickness', fontsize=13, fontweight='bold')
    ax.grid(True, alpha=0.3)

    # Panel C: f_cortex vs optical property variation (swept curves)
    ax = axes[1, 0]
    var_pct = opt_results['variations_pct']
    ax.plot(var_pct, opt_results['mua']['f_cortex'], 'b-o', linewidth=2,
            markersize=5, label='μa variation')
    ax.plot(var_pct, opt_results['musp']['f_cortex'], 'r-s', linewidth=2,
            markersize=5, label="μs' variation")
    ax.axhline(opt_results['baseline']['f_cortex'], color='k', linestyle='--',
               linewidth=1.5, alpha=0.5, label='Baseline')
    ax.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.4)
    ax.set_xlim([-32, 32])
    ax.set_xlabel('Optical Property Variation (%)', fontsize=12)
    ax.set_ylabel('Cortical Fraction (f_cortex)', fontsize=12)
    ax.set_title('(C) Sensitivity to Optical Property Variation', fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    # Panel D: κ(PV) vs optical property variation (swept curves)
    ax = axes[1, 1]
    ax.plot(var_pct, opt_results['mua']['kappa_pv'], 'b-o', linewidth=2,
            markersize=5, label='μa variation')
    ax.plot(var_pct, opt_results['musp']['kappa_pv'], 'r-s', linewidth=2,
            markersize=5, label="μs' variation")
    ax.axhline(opt_results['baseline']['kappa_pv'], color='k', linestyle='--',
               linewidth=1.5, alpha=0.5, label='Baseline')
    ax.axvline(0, color='gray', linestyle=':', linewidth=1, alpha=0.4)
    ax.set_xlim([-32, 32])
    ax.set_xlabel('Optical Property Variation (%)', fontsize=12)
    ax.set_ylabel('Partial Volume Correction (κ(PV))', fontsize=12)
    ax.set_title('(D) κ_PV Sensitivity to Optical Properties', fontsize=13, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"  Figure saved to: {save_path}")
    plt.close()
    return fig


# =============================================================================
# EXPANDED ANALYSIS FUNCTIONS
# =============================================================================

def analyze_HbR_performance(df: pd.DataFrame) -> Dict:
    """Analyze HbR correction performance at each SDS and overall."""
    print("\n" + "="*70)
    print("ANALYSIS A: HbR PERFORMANCE")
    print("="*70)

    results = {}
    rmse_mbll_hbr = np.sqrt((df['error_mbll_HbR']**2).mean())
    rmse_corr_hbr = np.sqrt((df['error_corr_HbR']**2).mean())
    improvement_hbr = (1 - rmse_corr_hbr/rmse_mbll_hbr) * 100

    print("\nOverall HbR RMSE:")
    print(f"  MBLL:      {rmse_mbll_hbr:.4f} µM")
    print(f"  Corrected: {rmse_corr_hbr:.4f} µM")
    print(f"  Improvement: {improvement_hbr:.1f}%")

    results['overall'] = {
        'rmse_mbll': rmse_mbll_hbr,
        'rmse_corrected': rmse_corr_hbr,
        'improvement_pct': improvement_hbr
    }

    print(f"\n{'SDS (mm)':<12} {'MBLL RMSE':<15} {'Corr RMSE':<15} {'Improv %':<12}")
    print("-" * 54)

    for sds in sorted(df['sds_mm'].unique()):
        subset = df[df['sds_mm'] == sds]
        rmse_mbll = np.sqrt((subset['error_mbll_HbR']**2).mean())
        rmse_corr = np.sqrt((subset['error_corr_HbR']**2).mean())
        improvement = (1 - rmse_corr/rmse_mbll) * 100

        print(f"{sds:<12.1f} {rmse_mbll:<15.4f} {rmse_corr:<15.4f} {improvement:<12.1f}")
        results[f'sds_{sds}'] = {
            'rmse_mbll': rmse_mbll,
            'rmse_corrected': rmse_corr, 'improvement_pct': improvement
        }

    return results


def analyze_per_subject_variability(df: pd.DataFrame) -> Dict:
    """Analyze per-subject variability of corrected estimates."""
    print("\n" + "="*70)
    print("ANALYSIS B: PER-SUBJECT VARIABILITY")
    print("="*70)

    results = {}
    print(f"\n{'SDS (mm)':<12} {'Mean β':<12} {'Std β':<12} {'Min Err':<12} {'Max Err':<12} {'Range':<12}")
    print("-" * 72)

    for sds in sorted(df['sds_mm'].unique()):
        subset = df[df['sds_mm'] == sds]
        mean_corr = subset['corrected_beta_HbO'].mean()
        std_corr = subset['corrected_beta_HbO'].std()
        min_err = subset['error_corr_HbO'].min()
        max_err = subset['error_corr_HbO'].max()
        err_range = max_err - min_err

        print(f"{sds:<12.1f} {mean_corr:<12.4f} {std_corr:<12.4f} {min_err:<12.4f} {max_err:<12.4f} {err_range:<12.4f}")
        results[f'sds_{sds}'] = {
            'mean_beta': mean_corr, 'std_beta': std_corr,
            'min_error': min_err, 'max_error': max_err, 'error_range': err_range
        }

    subject_mae = df.groupby('subject')['error_corr_HbO'].apply(lambda x: np.mean(np.abs(x)))
    print(f"\nPer-subject MAE: {subject_mae.mean():.4f} ± {subject_mae.std():.4f} µM")
    print(f"Range: {subject_mae.min():.4f} to {subject_mae.max():.4f} µM")

    results['overall_mae'] = {
        'mean': subject_mae.mean(), 'std': subject_mae.std(),
        'min': subject_mae.min(), 'max': subject_mae.max()
    }
    return results


def analyze_cortical_snr(data: Dict, wavelength_nm: int = 760) -> Dict:
    """Cortical optical-density SNR versus SDS (reproduces the SNR table).

    The cortical HbO2 contribution to the long-channel OD at a given wavelength
    is, in decadic (base-10) OD units,

        dOD_cortex(t) = eps_HbO2 * dHbO2_ctx(t) * L * DPF(lambda) * f_cortex(lambda),

    the f_cortex-weighted cortical HbO2 part of the forward signal model (Eq.
    signal_model), expressed in the same base-10 units as the assumed OD-noise level
    sigma_OD = 1e-3 (equivalently sigma_OD = ln(10) * 1e-3 in Napierian OD; the
    ln(10) cancels between signal and noise so the SNR is convention-independent).
    The HbO2 term alone is used because the short-SDS limit analysed in the text
    concerns HbO2 recovery; at 760 nm the HbR term partially cancels HbO2 at the
    peak, so the combined bracket would understate the HbO2 signal of interest.
    We report the peak |dOD_cortex| over the block-design time course, so
    SNR = peak(dOD_cortex) / sigma_OD.  This quantifies why kappa_PV = 1/f_cortex
    amplifies noise rather than signal at short SDS (SNR < 1), and becomes
    useful only where SNR >~ 2 (long SDS).  Values are averaged over subjects.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS: CORTICAL OD SNR VERSUS SEPARATION")
    print("=" * 70)
    extinction = ExtinctionCoefficients()
    E = extinction.get_matrix([wavelength_nm])[0]      # [eps_HbO2, eps_HbR] mM^-1 cm^-1 (decadic)
    eps_hbo2, eps_hbr = float(E[0]), float(E[1])
    dpf = scholkmann_dpf(wavelength_nm)
    sigma_OD = 1e-3                                     # assumed additive OD measurement-noise floor (fixed, separation-independent)

    print(f"\n  lambda = {wavelength_nm} nm | DPF = {dpf:.3f} | "
          f"sigma_OD = {sigma_OD:.3e} (base-10 OD)")
    print(f"\n{'SDS (mm)':<10}{'f_cortex':<12}{'dOD_cortex':<14}{'sigma_OD':<12}"
          f"{'SNR':<8}{'kappa_PV':<10}")
    print("-" * 66)

    results = {}
    # group SDS by channel; cortical ground truth is per subject (scaled), average it
    sds_values = sorted({ch['sds_mm'] for s in data['subjects'] for ch in s['channels']})
    for sds in sds_values:
        L_cm = sds / 10.0                              # mm -> cm
        peaks, fcs = [], []
        for subj in data['subjects']:
            for ch in subj['channels']:
                if ch['sds_mm'] != sds:
                    continue
                f_c = ch['f_cortex_dict'][wavelength_nm]
                # ground-truth cortical HbO2 is stored in micromolar; convert to
                # mM (1e-3) so eps [mM^-1 cm^-1] * C [mM] * L [cm] is a
                # dimensionless base-10 OD (HbO2 term only; see docstring).
                dod = (eps_hbo2 * ch['truth_HbO_cortex'] * 1e-3
                       * L_cm * dpf * f_c)
                peaks.append(float(np.max(np.abs(dod))))
                fcs.append(f_c)
        dod_peak = float(np.mean(peaks))
        f_c = float(np.mean(fcs))
        snr = dod_peak / sigma_OD
        kpv = 1.0 / f_c
        print(f"{sds:<10.0f}{f_c:<12.4f}{dod_peak:<14.3e}{sigma_OD:<12.3e}"
              f"{snr:<8.2f}{kpv:<10.1f}")
        results[f'sds_{sds}'] = dict(f_cortex=f_c, dOD_cortex=dod_peak,
                                     sigma_OD=sigma_OD, SNR=snr, kappa_PV=kpv)
    print("\n  Interpretation: correction helps only where SNR >~ 2 (long SDS);")
    print("  at 25 mm the cortical OD sits at/below the noise floor (SNR < 1).")
    return results


def analyze_thickness_robustness(sds_mm: float = 30.0, wavelength_nm: int = 760) -> Dict:
    """Analyze sensitivity to superficial layer thickness."""
    print("\n" + "="*70)
    print("ANALYSIS C: ROBUSTNESS TO SUPERFICIAL LAYER THICKNESS")
    print("="*70)

    thicknesses = [8, 9, 10, 11, 12, 13, 14, 15, 16]
    results = {}

    # f_cortex vs superficial thickness at SDS = 30 mm. PREFERRED SOURCE: the
    # versioned robustness_secondary.json produced by mc_robustness_sweeps.py, read
    # directly so the figure/table are an exact function of that artifact. If the
    # file is absent we fall back to the archived MC arrays (clearly labelled).
    _rob = _robustness_from_json()   # raises in RELEASE mode if missing/unverified
    if _rob is not None:
        _MC_T, _MC_F = _rob['thickness_T'], _rob['thickness_F']
        print(f"  (thickness sweep read from {_rob['source']})")
    else:
        import os as _os
        if _os.environ.get("FNIRS_ALLOW_FALLBACK") != "1":
            raise FileNotFoundError(
                "robustness_secondary.json missing/unverified. Run mc_robustness_sweeps.py, "
                "or set FNIRS_ALLOW_FALLBACK=1 to use the archived development arrays.")
        _MC_T = [8.0, 10.0, 12.0, 14.0, 16.0]
        _MC_F = [0.2392, 0.1217, 0.0626, 0.0301, 0.0151]   # archived DEV fallback only
        print("  [DEV] robustness_secondary.json not found; using archived MC arrays "
              "(FNIRS_ALLOW_FALLBACK=1)")
    def _f_mc(T):
        return float(np.exp(np.interp(float(T), _MC_T, np.log(_MC_F))))

    print(f"\n{'Thickness':<12} {'f_cortex':<15} {'κ(PV)':<15} {'Δf %':<15} {'Δκ %':<15}")
    print("-" * 72)

    f_baseline = None
    k_baseline = None

    for i, thickness in enumerate(thicknesses):
        f_cortex = _f_mc(thickness)
        kappa_pv = 1.0 / f_cortex

        if i == 0:
            f_baseline = f_cortex
            k_baseline = kappa_pv
            delta_f = 0.0
            delta_k = 0.0
        else:
            delta_f = (f_cortex - f_baseline) / f_baseline * 100
            delta_k = (kappa_pv - k_baseline) / k_baseline * 100

        print(f"{thickness:<12} {f_cortex:<15.4f} {kappa_pv:<15.4f} {delta_f:<15.2f} {delta_k:<15.2f}")
        results[thickness] = {
            'f_cortex': f_cortex, 'kappa_pv': kappa_pv,
            'delta_f_pct': delta_f, 'delta_kappa_pct': delta_k
        }

    return results


# Converged Monte-Carlo optical-property sensitivity of f_cortex / kappa_PV at
# 30 mm (code/mc_2layer.py, Henyey-Greenstein g=0.9, 3e5 photons, seed 1; the
# cortical-mu_a curve uses correlated reweighting of a single photon ensemble,
# so it is low-variance).  These REPLACE the earlier analytical-kernel sweep,
# which (i) used the non-converged Cartesian-trapezoid Hankel f_cortex (~0.29 at
# 30 mm, an order of magnitude too high; see fcortex_source / mc_production and Hankel
# convergence) and (ii) had its kappa column decoupled from f.  Grid is
# -30,-15,0,+15,+30 % perturbation of the cortical optical property.
# ARCHIVED MC results (two-layer white MC at SDS=30 mm, 760 nm): mu_a by correlated
# reweighting of one photon ensemble, mu_s' by independent runs.  Regenerated from
# first principles by mc_robustness_sweeps.py (-> results/robustness_secondary.json)
# with full provenance; kept here for the figure and consistent within MC noise.
_OPT_SENS_MC = {
    'variations_pct': [-30, -15, 0, 15, 30],
    # cortical mu_a perturbation (correlated estimator -> precise)
    'mua':  {'f': [0.0930, 0.0758, 0.0629, 0.0530, 0.0453],
             'k': [10.748, 13.187, 15.907, 18.883, 22.091]},
    # cortical mu_s' perturbation (independent runs -> ~few-% MC noise)
    'musp': {'f': [0.0712, 0.0650, 0.0629, 0.0582, 0.0525],
             'k': [14.037, 15.378, 15.907, 17.176, 19.030]},
    'baseline': {'f': 0.0629, 'k': 15.907},
}


def analyze_optical_property_robustness(sds_mm: float = 30.0, wavelength_nm: int = 760) -> Dict:
    """Sensitivity of f_cortex / kappa_PV to +/-30% variation in the cortical
    optical properties, from the converged two-layer Monte Carlo (mc_2layer.py)
    at 30 mm.  Reported as self-consistent f_cortex and kappa_PV = 1/f_cortex
    (the previous version mixed a non-converged analytical f_cortex with an
    MC kappa and is no longer used)."""
    print("\n" + "="*70)
    print("ANALYSIS D: ROBUSTNESS TO OPTICAL PROPERTY VARIATION (Monte Carlo)")
    print("="*70)

    # PREFERRED SOURCE: robustness_secondary.json (mc_robustness_sweeps.py), hash-
    # verified; RELEASE mode requires it. The archived arrays are a DEV-only fallback
    # (FNIRS_ALLOW_FALLBACK=1).
    _rob = _robustness_from_json()   # raises in RELEASE mode if missing/unverified
    if _rob is None:
        import os as _os
        if _os.environ.get("FNIRS_ALLOW_FALLBACK") != "1":
            raise FileNotFoundError(
                "robustness_secondary.json missing/unverified. Run mc_robustness_sweeps.py, "
                "or set FNIRS_ALLOW_FALLBACK=1 to use the archived development arrays.")
        print("  [DEV] using archived optical-sensitivity arrays (FNIRS_ALLOW_FALLBACK=1)")
    tbl = _rob['opt'] if _rob is not None else _OPT_SENS_MC
    if _rob is not None:
        print(f"  (optical sweep read from {_rob['source']})")
    variations_pct = np.array(tbl['variations_pct'], dtype=float)
    f_base = tbl['baseline']['f']
    k_base = tbl['baseline']['k']

    baseline_model = TwoLayerModel.default_adult(z_superficial_mm=12.0)
    props_base = baseline_model.props_cortex[wavelength_nm]
    print(f"\nBaseline (Monte Carlo, 30 mm): μa={props_base.mua:.4f}, "
          f"μs'={props_base.musp:.4f}, f_cortex={f_base:.4f}, κ_PV={k_base:.3f}")

    results: Dict = {'variations_pct': variations_pct}
    for prop, label in (('mua', 'μa'), ('musp', "μs'")):
        f_arr = np.array(tbl[prop]['f'], dtype=float)
        k_arr = np.array(tbl[prop]['k'], dtype=float)
        print(f"\n{label} Variation:")
        print(f"{'Variation':<12} {'f_cortex':<12} {'Δf %':<12} {'κ(PV)':<12} {'Δκ %':<12}")
        print("-" * 60)
        for v, f, k in zip(variations_pct, f_arr, k_arr):
            df = (f - f_base) / f_base * 100
            dk = (k - k_base) / k_base * 100
            print(f"{v:+.0f}%{'':<8} {f:<12.4f} {df:<+12.2f} {k:<12.3f} {dk:<+12.2f}")
        results[prop] = {
            'f_cortex': f_arr, 'kappa_pv': k_arr,
            'f_cortex_range': (float(f_arr.min()), float(f_arr.max())),
            'kappa_range': (float(k_arr.min()), float(k_arr.max())),
            'kappa_delta_pct': ((k_arr.min() - k_base) / k_base * 100,
                                (k_arr.max() - k_base) / k_base * 100),
        }
    results['baseline'] = {'f_cortex': f_base, 'kappa_pv': k_base}

    dmua = results['mua']['kappa_delta_pct']
    dmusp = results['musp']['kappa_delta_pct']
    print(f"\n  κ_PV change for ±30% cortical μa:  {dmua[0]:+.1f}% to {dmua[1]:+.1f}%")
    print(f"  κ_PV change for ±30% cortical μs': {dmusp[0]:+.1f}% to {dmusp[1]:+.1f}%")
    print("  (Monte Carlo; μa dominates. Still smaller than the >100% change "
          "from a ±2 mm thickness error.)")
    return results


def analyze_grid_convergence(sds_mm: float = 30.0, wavelength_nm: int = 760) -> Dict:
    """
    Grid convergence study for numerical stability verification.

    Varies dx (x-resolution) and dz (depth resolution) on a 4×4 grid,
    with the y-dimension fixed at its default resolution (dy_mm=2.0 mm,
    y_range_mm=50.0 mm) since the 3D integral's y-axis is separately
    verified by the y-symmetry construction.

    The finest (dx=0.25, dz=0.1) grid is used as the reference solution.
    """
    print("\n" + "="*70)
    print("ANALYSIS E: GRID CONVERGENCE")
    print("="*70)

    model = TwoLayerModel.default_adult(z_superficial_mm=12.0)
    results = {}

    # Fixed y-grid parameters (not varied in this convergence study)
    y_range_fixed = 50.0
    dy_fixed = 2.0
    n_y_fixed = int(y_range_fixed / dy_fixed) + 1  # 26 unique y-slices (with symmetry)

    dx_values = [2.0, 1.0, 0.5, 0.25]
    dz_values = [1.0, 0.5, 0.25, 0.1]

    # First pass: compute f_cortex for all grid combinations
    grid_results = {}
    for dx in dx_values:
        for dz in dz_values:
            _, _, f = compute_sensitivity_map(
                sds_mm, wavelength_nm, model,
                dz_mm=dz, dx_mm=dx,
                y_range_mm=y_range_fixed, dy_mm=dy_fixed
            )
            n_z = int(35.0 / dz)
            n_x = int(100.0 / dx) + 1
            n_pts = n_z * n_y_fixed * n_x
            grid_results[(dx, dz)] = {'f_cortex': f, 'n_points': n_pts}

    # Use the finest grid (smallest dx, smallest dz) as reference
    f_finest = grid_results[(min(dx_values), min(dz_values))]['f_cortex']

    print(f"\nFixed y-grid: dy={dy_fixed} mm, y_range={y_range_fixed} mm ({n_y_fixed} slices + symmetry)")
    print(f"\n{'dx (mm)':<12} {'dz (mm)':<12} {'f_cortex':<15} {'Δf from finest':<16} {'N_points':<12}")
    print("-" * 67)

    for dx in dx_values:
        for dz in dz_values:
            f = grid_results[(dx, dz)]['f_cortex']
            n_pts = grid_results[(dx, dz)]['n_points']
            delta = (f - f_finest) / f_finest * 100 if f_finest > 0 else 0.0

            print(f"{dx:<12.2f} {dz:<12.2f} {f:<15.6f} {delta:<16.4f} {n_pts:<12}")
            results[f'dx{dx}_dz{dz}'] = {'f_cortex': f, 'delta_f_pct': delta, 'n_points': n_pts}

    return results


def benchmark_computation_time(sds_mm: float = 30.0, wavelength_nm: int = 760,
                               n_trials: int = 3) -> Dict:
    """Benchmark single-channel sensitivity computation time."""
    print("\n" + "="*70)
    print("ANALYSIS F: COMPUTATION TIME BENCHMARKING")
    print("="*70)

    model = TwoLayerModel.default_adult(z_superficial_mm=12.0)
    times = []

    print("\n  Benchmarking 3D sensitivity computation (dx=1.0, dz=0.5, dy=2.0 mm)...")
    for trial in range(n_trials):
        start = time.time()
        compute_sensitivity_map(sds_mm, wavelength_nm, model, dx_mm=1.0, dz_mm=0.5)
        elapsed = time.time() - start
        times.append(elapsed)
        print(f"    Trial {trial+1}: {elapsed:.4f} s")

    times = np.array(times)
    print(f"\n  Mean: {times.mean():.4f} ± {times.std():.4f} s")

    n_channels = 40
    n_wavelengths = 2  # two measurement wavelengths (760 and 850 nm)
    total_est = times.mean() * n_channels * n_wavelengths
    print(f"  Est. full pipeline ({n_channels} ch x {n_wavelengths} wl): {total_est:.2f} s")

    return {
        'mean_time': float(times.mean()),
        'std_time': float(times.std()),
        'n_channels': n_channels,
        'n_wavelengths': n_wavelengths,
        'total_estimate': float(total_est)
    }


# =============================================================================
# FINITE-DIFFERENCE AND MULTI-LAYER ANALYSES (NEW)
# =============================================================================

# NOTE: An earlier version of this file contained a function named
# ``solve_diffusion_fd`` that was presented as an independent finite-difference
# cross-check of f_cortex.  It did not actually assemble or solve a
# finite-difference system: it filled the source/detector fields with
# homogeneous-medium Green's functions exp(-mu_eff*r)/(4*pi*D*r) and integrated
# their product.  Because it was not a genuine FD solve and its results were not
# reproduced by any runnable command, that function and the associated
# "finite-difference / adaptive-Kienle validation" table have been REMOVED.
# The converged f_cortex is validated by (a) Monte-Carlo convergence
# (L_max / z_max / batch checks in mc_production.py) and (b) agreement with
# published atlas Monte Carlo (Strangman, Li & Zhang, Colin27) -- see the
# manuscript's atlas-comparison paragraph.  No independent FD/adaptive-Kienle
# numerical validation is claimed.


# ============================================================================
# CSF (THREE-LAYER) CORTICAL SENSITIVITY  --  MONTE-CARLO CALIBRATED
# ============================================================================
# An earlier diffusion finite-difference approximation estimated the effect of an
# explicit cerebrospinal-fluid (CSF) layer and returned a ~99% COLLAPSE of cortical
# sensitivity.  That result is a numerical artifact: the diffusion approximation is
# invalid in the thin (1-2 mm) CSF layer because that layer is much thinner than its
# transport mean free path (1/mu_s' ~ 4 mm at mu_s' ~ 0.25 mm^-1), so photons cross
# it in far less than one transport length and undergo too few direction-randomizing
# events for diffusion theory to hold; a 99% collapse is also physically implausible
# (fNIRS detects cortical activation through real CSF).
#
# It is replaced by the anisotropic (Henyey-Greenstein g=0.9) white Monte Carlo of
# the PRODUCTION source (code/mc_production.py, read via fcortex_source.py), which
# shows that inserting a low-scattering CSF layer INCREASES the cortical sensitivity
# fraction (the CSF 'light-piping' effect), by a factor
#     gamma(SDS) = f_cortex(3-layer) / f_cortex(2-layer)
# of ~1.4-1.8 (decreasing with SDS), in agreement with the light-transport literature.
# The gamma values come from fcortex_source.gamma_csf (production Monte Carlo); the
# legacy isotropic mc_csf.py is NOT used here.
def csf_gamma(sds_mm: float) -> float:
    """CSF amplification f_cortex(3L)/f_cortex(2L) from the single production
    Monte-Carlo source (760 nm; used only for the illustrative CSF/mismatch
    analysis).  No hard-coded table."""
    return _fcs.gamma_csf(sds_mm, 760)

def f_cortex_three_layer(sds_mm: float, f_cortex_2l: float) -> float:
    """Three-layer (with CSF) cortical sensitivity fraction, Monte-Carlo calibrated."""
    return float(min(csf_gamma(sds_mm) * f_cortex_2l, 0.999))


def analyze_csf_layer_effect(sds_list: List[float] = [25, 35]) -> Dict:
    """
    Analyze cortical sensitivity with explicit CSF layer (three-layer model).

    Three-layer model:
      Layer 1 (scalp-skull): 10 mm, μ_a=0.015, μ_s'=1.0
      Layer 2 (CSF):         2 mm, μ_a=0.004, μ_s'=0.25
      Layer 3 (cortex):      semi-infinite, μ_a=0.020, μ_s'=0.80

    Compares against two-layer baseline (12mm superficial, cortex).
    """
    print("\n" + "="*70)
    print("ANALYSIS G: THREE-LAYER MODEL AND CSF EFFECTS")
    print("="*70)

    wavelength_nm = 760
    results = []

    # Two-layer baseline: use the converged Monte-Carlo f_cortex (matching
    # Table tab:results), NOT the non-converged analytical Hankel sensitivity,
    # so the absolute two-layer fractions printed here are consistent with the
    # rest of the manuscript.  The CSF effect is the Monte-Carlo ratio gamma.
    print("\n  Two-layer baseline (Monte-Carlo f_cortex; mc_2layer.py)...")

    for sds in sds_list:
        f_cortex_2l = f_cortex_mc(sds, wavelength_nm)
        print(f"    SDS={sds} mm: f_cortex(2L) = {f_cortex_2l:.4f}")

        # Three-layer (with CSF): production Monte-Carlo-calibrated amplification
        # (fcortex_source.gamma_csf, from mc_production.py; not the legacy mc_csf.py)
        f_cortex_3l = f_cortex_three_layer(sds, f_cortex_2l)

        increase_pct = (f_cortex_3l / f_cortex_2l - 1.0) * 100.0 if f_cortex_2l > 0 else 0.0

        print(f"    SDS={sds} mm: f_cortex(3L) = {f_cortex_3l:.4f}, increase = {increase_pct:.1f}%")

        results.append({
            'sds_mm': sds,
            'f_cortex_2l': f_cortex_2l,
            'f_cortex_3l': f_cortex_3l,
            'increase_pct': increase_pct
        })

    print("\n  Three-layer CSF Effect Summary:")
    print("  " + "-"*60)
    print(f"  {'SDS (mm)':>12} {'f_cortex(2L)':>15} {'f_cortex(3L)':>15} {'Increase (%)':>15}")
    print("  " + "-"*60)
    for r in results:
        print(f"  {r['sds_mm']:>12.0f} {r['f_cortex_2l']:>15.4f} {r['f_cortex_3l']:>15.4f} {r['increase_pct']:>15.1f}")
    print("  " + "-"*60)

    return {'results': results}


def analyze_model_mismatch(sds_mm: float = 30.0,
                           wavelengths: List[int] = [760, 850]) -> Dict:
    """
    Test correction performance under model mismatch using the full
    two-wavelength MBLL inversion (consistent with the main pipeline).

    Matched:    two-layer data corrected with the two-layer kappa_PV(lambda).
    Mismatched: three-layer (CSF) data corrected with the two-layer kappa_PV(lambda).

    Both HbO2 and HbR are recovered by a proper 2x2 MBLL inversion of the
    PV-corrected per-wavelength Delta OD; the partial-volume correction is the
    per-wavelength reciprocal cortical sensitivity 1/f_cortex(lambda).  This
    isolates the effect of using an incorrect tissue geometry (f_cortex) on the
    quantitative concentration estimate.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS H: MODEL MISMATCH TEST (2L vs 3L)")
    print("=" * 70)

    print("\n  Generating synthetic fNIRS data (n_samples=500)...")
    np.random.seed(42)
    n_samples = 500

    # Task regressor: two blocks of activity
    task_regressor = np.zeros(n_samples)
    task_regressor[100:200] = 1.0
    task_regressor[300:400] = 1.0
    task_regressor = (task_regressor - task_regressor.mean()) / (task_regressor.std() + 1e-10)

    # Ground-truth cortical concentrations (uM)
    hbo_true = 2.5 * task_regressor + 0.1 * np.random.randn(n_samples)
    hbr_true = -0.8 * task_regressor + 0.05 * np.random.randn(n_samples)

    extinction = ExtinctionCoefficients()
    dpf_dict = get_dpf_dict(wavelengths)

    # Wavelength-specific cortical sensitivity: two-layer (analytic Kienle)
    f2l = {wl: f_cortex_mc(sds_mm, wl) for wl in wavelengths}  # stable MC f_cortex

    # Three-layer (with CSF) cortical sensitivity, Monte-Carlo calibrated
    # (from the production CSF ratio gamma via fcortex_source):
    # CSF light-piping INCREASES f_cortex, so f3l > f2l and a two-layer model
    # UNDERESTIMATES the true cortical sensitivity (-> the two-layer kappa_PV overcorrects).
    f3l = {wl: f_cortex_three_layer(sds_mm, f2l[wl]) for wl in wavelengths}
    f3l_val = float(np.mean(list(f3l.values())))

    print("  f_cortex(2L): " + ", ".join(f"{wl}nm={f2l[wl]:.4f}" for wl in wavelengths))
    print(f"  f_cortex(3L): {f3l_val:.4f} (Monte-Carlo calibrated; CSF light-piping)")

    sigma_od = 0.001  # measurement noise std at OD level

    def _forward(fdict):
        """Generate per-wavelength Delta OD matrix (n_samples, n_wl)."""
        od = np.zeros((n_samples, len(wavelengths)))
        for i, wl in enumerate(wavelengths):
            eps_hbo = extinction.coefficients[wl][0] * 1e-4   # mm^-1 / uM
            eps_hbr = extinction.coefficients[wl][1] * 1e-4
            clean = (eps_hbo * hbo_true + eps_hbr * hbr_true) * sds_mm * dpf_dict[wl] * fdict[wl]
            od[:, i] = clean + sigma_od * np.random.randn(n_samples)
        return od

    def _invert(od_matrix, f_corr=None):
        """Two-wavelength MBLL inversion; f_corr applies OD-level PV correction."""
        return mbll_inversion_batch(od_matrix, wavelengths, sds_mm, dpf_dict,
                                    extinction, f_cortex_dict=f_corr)

    def _rmse(est, truth):
        return float(np.sqrt(np.mean((est - truth) ** 2)))

    # --- Scenario 1: matched (2L data -> 2L correction) ---
    print("\n  Scenario 1: Matched (2L data -> 2L correction)")
    od_2l = _forward(f2l)
    mbll_hbo_2l, mbll_hbr_2l = _invert(od_2l, None)
    corr_hbo_2l, corr_hbr_2l = _invert(od_2l, f2l)

    rmse_mbll_hbo_2l = _rmse(mbll_hbo_2l, hbo_true)
    rmse_corr_hbo_2l = _rmse(corr_hbo_2l, hbo_true)
    improvement_hbo_2l = (1 - rmse_corr_hbo_2l / rmse_mbll_hbo_2l) * 100.0
    rmse_mbll_hbr_2l = _rmse(mbll_hbr_2l, hbr_true)
    rmse_corr_hbr_2l = _rmse(corr_hbr_2l, hbr_true)
    improvement_hbr_2l = (1 - rmse_corr_hbr_2l / rmse_mbll_hbr_2l) * 100.0
    print(f"    HbO2 RMSE: {rmse_mbll_hbo_2l:.3f} -> {rmse_corr_hbo_2l:.3f} ({improvement_hbo_2l:.1f}% improvement)")
    print(f"    HbR  RMSE: {rmse_mbll_hbr_2l:.3f} -> {rmse_corr_hbr_2l:.3f} ({improvement_hbr_2l:.1f}% improvement)")

    # --- Scenario 2: mismatched (3L data -> 2L correction) ---
    print("\n  Scenario 2: Mismatched (3L data -> 2L correction)")
    od_3l = _forward(f3l)
    mbll_hbo_3l, mbll_hbr_3l = _invert(od_3l, None)
    corr_hbo_mis, corr_hbr_mis = _invert(od_3l, f2l)

    rmse_mbll_hbo_3l = _rmse(mbll_hbo_3l, hbo_true)
    rmse_corr_hbo_mis = _rmse(corr_hbo_mis, hbo_true)
    improvement_hbo_mis = (1 - rmse_corr_hbo_mis / rmse_mbll_hbo_3l) * 100.0
    rmse_mbll_hbr_3l = _rmse(mbll_hbr_3l, hbr_true)
    rmse_corr_hbr_mis = _rmse(corr_hbr_mis, hbr_true)
    improvement_hbr_mis = (1 - rmse_corr_hbr_mis / rmse_mbll_hbr_3l) * 100.0
    print(f"    HbO2 RMSE: {rmse_mbll_hbo_3l:.3f} -> {rmse_corr_hbo_mis:.3f} ({improvement_hbo_mis:.1f}% {'improvement' if improvement_hbo_mis > 0 else 'worsening'})")
    print(f"    HbR  RMSE: {rmse_mbll_hbr_3l:.3f} -> {rmse_corr_hbr_mis:.3f} ({improvement_hbr_mis:.1f}% {'improvement' if improvement_hbr_mis > 0 else 'worsening'})")

    print("\n  Model Mismatch Summary:")
    print("  " + "-" * 70)
    print(f"  {'Scenario':<20} {'Metric':<12} {'MBLL RMSE':<15} {'Corrected RMSE':<15} {'Improvement':<12}")
    print("  " + "-" * 70)
    print(f"  {'Matched (2L->2L)':<20} {'HbO2':<12} {rmse_mbll_hbo_2l:>14.3f} {rmse_corr_hbo_2l:>14.3f} {improvement_hbo_2l:>11.1f}%")
    print(f"  {'':<20} {'HbR':<12} {rmse_mbll_hbr_2l:>14.3f} {rmse_corr_hbr_2l:>14.3f} {improvement_hbr_2l:>11.1f}%")
    print(f"  {'Mismatched (3L->2L)':<20} {'HbO2':<12} {rmse_mbll_hbo_3l:>14.3f} {rmse_corr_hbo_mis:>14.3f} {improvement_hbo_mis:>11.1f}%")
    print(f"  {'':<20} {'HbR':<12} {rmse_mbll_hbr_3l:>14.3f} {rmse_corr_hbr_mis:>14.3f} {improvement_hbr_mis:>11.1f}%")
    print("  " + "-" * 70)

    return {
        'f_cortex_2l': f2l, 'f_cortex_3l': f3l_val,
        'matched': {
            'rmse_mbll_hbo': rmse_mbll_hbo_2l, 'rmse_corr_hbo': rmse_corr_hbo_2l,
            'improvement_hbo': improvement_hbo_2l,
            'rmse_mbll_hbr': rmse_mbll_hbr_2l, 'rmse_corr_hbr': rmse_corr_hbr_2l,
            'improvement_hbr': improvement_hbr_2l,
        },
        'mismatched': {
            'rmse_mbll_hbo': rmse_mbll_hbo_3l, 'rmse_corr_hbo': rmse_corr_hbo_mis,
            'improvement_hbo': improvement_hbo_mis,
            'rmse_mbll_hbr': rmse_mbll_hbr_3l, 'rmse_corr_hbr': rmse_corr_hbr_mis,
            'improvement_hbr': improvement_hbr_mis,
        }
    }


def analyze_fcortex_sensitivity(sds_mm: float = 30.0,
                                wavelengths: List[int] = [760, 850],
                                perturbations: List[int] = [-30, -20, -10, 0, 10, 20, 30]) -> Dict:
    """
    Sensitivity analysis: how correction performance varies with f_cortex error,
    using the full two-wavelength MBLL inversion (matched two-layer scenario).

    The assumed cortical sensitivity is perturbed by a fixed fractional error
    applied to both wavelengths, and the resulting RMSE improvement relative to
    uncorrected MBLL is evaluated for HbO2 and HbR.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS I: f_cortex SENSITIVITY PERTURBATION")
    print("=" * 70)

    np.random.seed(42)
    n_samples = 500

    task_regressor = np.zeros(n_samples)
    task_regressor[100:200] = 1.0
    task_regressor[300:400] = 1.0
    task_regressor = (task_regressor - task_regressor.mean()) / (task_regressor.std() + 1e-10)

    hbo_true = 2.5 * task_regressor + 0.1 * np.random.randn(n_samples)
    hbr_true = -0.8 * task_regressor + 0.05 * np.random.randn(n_samples)

    extinction = ExtinctionCoefficients()
    dpf_dict = get_dpf_dict(wavelengths)

    # True per-wavelength cortical sensitivity
    f_true = {wl: f_cortex_mc(sds_mm, wl) for wl in wavelengths}  # stable MC f_cortex
    sigma_od = 0.001

    # Forward model with the true f_cortex
    od = np.zeros((n_samples, len(wavelengths)))
    for i, wl in enumerate(wavelengths):
        eps_hbo = extinction.coefficients[wl][0] * 1e-4
        eps_hbr = extinction.coefficients[wl][1] * 1e-4
        clean = (eps_hbo * hbo_true + eps_hbr * hbr_true) * sds_mm * dpf_dict[wl] * f_true[wl]
        od[:, i] = clean + sigma_od * np.random.randn(n_samples)

    def _rmse(est, truth):
        return float(np.sqrt(np.mean((est - truth) ** 2)))

    # Uncorrected MBLL (two-wavelength)
    mbll_hbo, mbll_hbr = mbll_inversion_batch(od, wavelengths, sds_mm, dpf_dict,
                                              extinction, f_cortex_dict=None)
    rmse_mbll_hbo = _rmse(mbll_hbo, hbo_true)
    rmse_mbll_hbr = _rmse(mbll_hbr, hbr_true)

    f_true_mean = float(np.mean(list(f_true.values())))
    print(f"\n  True f_cortex (mean) = {f_true_mean:.4f}")
    print(f"  Uncorrected MBLL RMSE: HbO2 = {rmse_mbll_hbo:.3f}, HbR = {rmse_mbll_hbr:.3f}")

    print("\n  Perturbing f_cortex assumption...")
    results = []
    print("\n  f_cortex Sensitivity Summary:")
    print("  " + "-" * 80)
    print(f"  {'Error (%)':>10} {'f_cortex_assumed':>18} {'HbO2 Improvement':>18} {'HbR Improvement':>18}")
    print("  " + "-" * 80)

    for pert in perturbations:
        f_assumed = {wl: f_true[wl] * (1.0 + pert / 100.0) for wl in wavelengths}
        corr_hbo, corr_hbr = mbll_inversion_batch(od, wavelengths, sds_mm, dpf_dict,
                                                  extinction, f_cortex_dict=f_assumed)
        rmse_corr_hbo = _rmse(corr_hbo, hbo_true)
        rmse_corr_hbr = _rmse(corr_hbr, hbr_true)
        improvement_hbo = (1 - rmse_corr_hbo / rmse_mbll_hbo) * 100.0
        improvement_hbr = (1 - rmse_corr_hbr / rmse_mbll_hbr) * 100.0
        f_assumed_mean = float(np.mean(list(f_assumed.values())))
        print(f"  {pert:>10d} {f_assumed_mean:>18.4f} {improvement_hbo:>18.1f}% {improvement_hbr:>18.1f}%")
        results.append({
            'error_pct': pert,
            'f_cortex_assumed': f_assumed_mean,
            'improvement_hbo': improvement_hbo,
            'improvement_hbr': improvement_hbr,
            'rmse_corr_hbo': rmse_corr_hbo,
            'rmse_corr_hbr': rmse_corr_hbr,
        })

    print("  " + "-" * 80)
    return {'results': results, 'f_cortex_true': f_true_mean}


def multiseed_operating_regime(n_seeds: int = 30, base_seed: int = 1000,
                               sds_list=(25.0, 30.0, 35.0, 38.0, 40.0)) -> List[Dict]:
    """Repeat the HbO2 RMSE-improvement evaluation over many independently
    GENERATED five-subject synthetic COHORTS.  NOTE: each seed regenerates the
    entire cohort -- subject amplitudes, systemic oscillations, phases,
    superficial signals AND the additive OD noise -- so this quantifies
    variability across independent synthetic-cohort realizations, not sensitivity
    to the measurement-noise draw alone.  The single-seed (seed 42) numbers are
    one such realization; this reports mean +/- SD and the fraction of cohorts
    improving at each separation, so the crossover can be described as conditional
    on the assumed signal-and-noise model rather than a hard threshold.

    Returns a list of dicts (one per SDS) with mean, sd and fraction-improving.
    """
    per_sds = {s: [] for s in sds_list}
    for k in range(n_seeds):
        data = generate_synthetic_fnirs_data(
            n_subjects=5, duration_s=660.0, sfreq=7.8125,
            sds_long_mm=list(sds_list), sds_short_mm=10.0,
            wavelengths=[760, 850], seed=base_seed + k)
        df = run_analysis_pipeline(data, n_s=500)
        for s in sds_list:
            sub = df[np.isclose(df['sds_mm'], s)]
            rm = np.sqrt((sub['error_mbll_HbO'] ** 2).mean())
            rc = np.sqrt((sub['error_corr_HbO'] ** 2).mean())
            per_sds[s].append((1.0 - rc / rm) * 100.0 if rm > 0 else np.nan)
    rows = []
    for s in sds_list:
        a = np.array(per_sds[s], float); a = a[~np.isnan(a)]
        rows.append(dict(sds_mm=float(s), n_seeds=int(a.size),
                         mean_improvement=float(a.mean()),
                         sd_improvement=float(a.std(ddof=1)) if a.size > 1 else 0.0,
                         frac_improving=float(np.mean(a > 0))))
    return rows


# =============================================================================
# MAIN FUNCTION
# =============================================================================

def main():
    """Main execution function demonstrating the complete pipeline."""
    print("="*70)
    print("fNIRS MBLL Bias Quantification and κ Correction Study")
    print("Two-layer Kienle SDA fluence | 3D sensitivity | all fixes applied")
    print("="*70)
    print()

    # =========================================================================
    # STEP 1: Generate synthetic data
    # =========================================================================
    print("Step 1: Generating synthetic fNIRS data...")
    print("-" * 50)

    data = generate_synthetic_fnirs_data(
        n_subjects=5,
        duration_s=660.0,
        sfreq=7.8125,
        sds_long_mm=[25.0, 30.0, 35.0, 38.0, 40.0],
        sds_short_mm=10.0,
        wavelengths=[760, 850],
        seed=42
    )

    print(f"  Generated data for {data['n_subjects']} subjects")
    print(f"  Duration: {data['n_samples']/data['sfreq']:.0f} s | Rate: {data['sfreq']} Hz")
    print(f"  Wavelengths: {data['wavelengths']} nm")
    print(f"  Long SDS: {data['sds_long_mm']} mm | Short SDS: {data['sds_short_mm']} mm")
    dpf_vals = get_dpf_dict(data['wavelengths'])
    for wl, dpf in dpf_vals.items():
        print(f"  DPF({wl} nm) = {dpf:.3f}  (Scholkmann & Wolf 2013)")
    print()

    # =========================================================================
    # STEP 2: Run analysis pipeline
    # =========================================================================
    print("Step 2: Running analysis pipeline (3D Kienle fluence + all fixes)...")
    print("-" * 50)

    results_df = run_analysis_pipeline(data, n_s=500)
    print(f"  Analyzed {len(results_df)} channel-subject combinations")
    print()

    # =========================================================================
    # STEP 3: Core HbO₂ results
    # =========================================================================
    print("Step 3: Core HbO₂ Results")
    print("-" * 50)

    print("\nCorrection Factors by SDS:")
    summary = results_df.groupby('sds_mm').agg({
        'f_cortex_computed': 'mean',
        'kappa_PV': 'mean',
        'V_SSR_760': 'mean',
        'V_SSR_850': 'mean',
        'error_mbll_HbO': lambda x: np.sqrt((x**2).mean()),
        'error_corr_HbO': lambda x: np.sqrt((x**2).mean()),
    }).rename(columns={
        'f_cortex_computed': 'f_cortex',
        'V_SSR_760': 'V_SSR_760',
        'V_SSR_850': 'V_SSR_850',
        'error_mbll_HbO': 'RMSE_MBLL',
        'error_corr_HbO': 'RMSE_Corrected'
    })
    print(summary.round(3).to_string())

    rmse_mbll = np.sqrt((results_df['error_mbll_HbO']**2).mean())
    rmse_corr = np.sqrt((results_df['error_corr_HbO']**2).mean())
    improvement = (1 - rmse_corr/rmse_mbll) * 100

    print("\nOverall HbO₂ Performance:")
    print(f"  MBLL RMSE:      {rmse_mbll:.4f} µM")
    print(f"  Corrected RMSE: {rmse_corr:.4f} µM")
    print(f"  Improvement:    {improvement:.1f}%")

    print("\nκ Factor Statistics (mean ± std):")
    print(f"  κ(DPF):   {results_df['kappa_DPF'].mean():.3f} ± {results_df['kappa_DPF'].std():.3f}")
    print(f"  κ(PV):    {results_df['kappa_PV'].mean():.3f} ± {results_df['kappa_PV'].std():.3f}  (applied)")
    # V_SSR is a per-wavelength diagnostic (NOT applied).  No cross-wavelength
    # mean or κ_total is reported, because the 760/850 nm values differ too much
    # for their arithmetic mean to be a meaningful quantity.
    print("  V(SSR) per wavelength [diagnostic, not applied]:")
    wl_cols = [c for c in results_df.columns if c.startswith('V_SSR_')]
    for col in sorted(wl_cols):
        wl = col.replace('V_SSR_', '')
        r2_col = f'R2_SS_{wl}'
        print(f"    λ={wl} nm: V_SSR = {results_df[col].mean():.3f} ± {results_df[col].std():.3f}"
              f"  [R²_SS = {results_df[r2_col].mean():.3f}]")
    print()

    # =========================================================================
    # STEP 4: Expanded analyses
    # =========================================================================
    print("Step 4: Expanded Analyses")
    print("-" * 50)

    hbr_results    = analyze_HbR_performance(results_df)
    subj_results   = analyze_per_subject_variability(results_df)
    snr_results    = analyze_cortical_snr(data, wavelength_nm=760)
    thickness_results = analyze_thickness_robustness()
    opt_results    = analyze_optical_property_robustness()
    analyze_grid_convergence()
    time_results   = benchmark_computation_time()

    # NOTE: the secondary robustness sweeps (f_cortex vs superficial thickness and
    # vs cortical mu_a/mu_s') are NO LONGER serialized here.  Stamping the embedded
    # archived arrays (_MC_F / _OPT_SENS_MC below) with this script's current run
    # metadata would misattribute provenance -- it would record when the constants
    # were written out, not when the Monte-Carlo sweeps that produced them ran.
    # Those sweeps are now regenerated FROM FIRST PRINCIPLES by a dedicated,
    # executable script (mc_robustness_sweeps.py), which writes
    # results/robustness_secondary.json with genuine per-configuration provenance
    # (seed, photon count, geometry, optical properties, git commit, command,
    # timestamp).  The plotted values used here (analyze_thickness_robustness /
    # analyze_optical_property_robustness) are archived MC results and agree with
    # that fresh run at the few-percent level expected between independent MC runs.
    print("  Secondary robustness sweeps: see mc_robustness_sweeps.py -> "
          "results/robustness_secondary.json (executable, first-principles MC).")

    # =========================================================================
    # STEP 4B: Multi-layer and validation analyses (NEW)
    # =========================================================================
    print("\nStep 4B: Multi-layer and Validation Analyses")
    print("-" * 50)

    analyze_csf_layer_effect()
    # Evaluate the model-mismatch and f_cortex-sensitivity analyses at the
    # long-separation operating point (38 mm; Section "Operating regime"), so
    # the console output matches the manuscript Tables for these analyses.
    # (Previously these used the 30 mm default, which did not reproduce the
    # reported 38 mm values.)
    analyze_model_mismatch(sds_mm=38.0)
    analyze_fcortex_sensitivity(sds_mm=38.0)

    # =========================================================================
    # STEP 5: Generate all figures
    # =========================================================================
    print("\n" + "="*70)
    print("Step 5: Generating figures...")
    print("-" * 50)

    plot_time_series_example(data, results_df, subject_idx=0, channel_idx=2,
                                     save_path='figure1_timeseries.png')
    plot_results(results_df, save_path='figure2_summary.png')
    plot_robustness(thickness_results, opt_results,
                           save_path='figure3_robustness.png')

    # =========================================================================
    # STEP 6: Save results
    # =========================================================================
    print("\nStep 6: Saving results...")
    print("-" * 50)

    results_df.to_csv('fnirs_kappa_results.csv', index=False)
    print("  Saved: fnirs_kappa_results.csv")

    summary.to_csv('fnirs_kappa_summary_by_sds.csv')
    print("  Saved: fnirs_kappa_summary_by_sds.csv")

    # --- Operating-regime variability across independent synthetic cohorts ---
    print("\n" + "="*70)
    print("ANALYSIS I: SYNTHETIC-COHORT OPERATING REGIME")
    print("  (HbO2 RMSE improvement across 30 independently generated 5-subject cohorts)")
    print("="*70)
    ms_rows = multiseed_operating_regime(n_seeds=30, base_seed=1000)
    print(f"\n{'SDS (mm)':<10}{'Mean impr.':<13}{'SD':<9}{'Frac. cohorts improving':<24}")
    print("-"*56)
    for r in ms_rows:
        print(f"{r['sds_mm']:<10.0f}{r['mean_improvement']:<+13.1f}"
              f"{r['sd_improvement']:<9.1f}{100*r['frac_improving']:<24.0f}")
    json.dump(dict(
        description="HbO2 RMSE improvement across 30 independently generated 5-subject "
                    "synthetic cohorts (each seed regenerates amplitudes, systemics, phases, "
                    "superficial signals and additive OD noise); not a noise-only analysis.",
        base_seed=1000, n_cohorts=30, per_sds=ms_rows,
        _meta=_synth_provenance("fnirs_kappa_synthetic_validation.py:multiseed_operating_regime")),
        open('multiseed_operating_regime.json', 'w'), indent=1)
    print("  Saved: multiseed_operating_regime.json")

    print()
    print("="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)
    print("\nKey Results:")
    print(f"  HbO₂ RMSE improvement: {improvement:.1f}%")
    print(f"  HbR  RMSE improvement: {hbr_results['overall']['improvement_pct']:.1f}%")
    print(f"  Per-subject MAE: {subj_results['overall_mae']['mean']:.4f} ± {subj_results['overall_mae']['std']:.4f} µM")
    print(f"  κ(PV) mean (applied):  {results_df['kappa_PV'].mean():.3f}")
    print(f"  V(SSR) [diagnostic, per λ]:  "
          f"760 nm {results_df['V_SSR_760'].mean():.3f}, "
          f"850 nm {results_df['V_SSR_850'].mean():.3f}")
    print(f"  Computation time: {time_results['mean_time']:.4f} s/channel")
    print()

    return data, results_df


if __name__ == "__main__":
    data, results = main()
