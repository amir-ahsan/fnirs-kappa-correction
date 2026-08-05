#!/usr/bin/env python3
"""
Single source of truth for cortical sensitivity fractions.

Both the synthetic-validation pipeline and the real-data pipeline import their
f_cortex values from here, which reads the versioned production Monte-Carlo file
``fcortex_production.json`` (written by mc_production.py). There are NO hard-coded
sensitivity tables and NO 760->850 ratio shortcut anywhere downstream: the
wavelength-specific two-layer fractions and the wavelength-specific CSF ratios
gamma(lambda, SDS) come directly from that one file.

If the production file is absent, a clear error is raised so that no analysis can
silently fall back to stale values.
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np

_PROD_NAME = "fcortex_production.json"


def _find_production() -> Path:
    here = Path(__file__).resolve().parent
    for cand in (here / _PROD_NAME, Path.cwd() / _PROD_NAME):
        if cand.exists():
            return cand
    raise FileNotFoundError(
        f"{_PROD_NAME} not found next to fcortex_source.py or in the CWD. "
        f"Generate it first:  python mc_production.py -N 2000000 --batches 16 "
        f"--out {Path(__file__).resolve().parent / 'fcortex_production'}")


_DATA = json.load(open(_find_production()))
WAVELENGTHS = sorted(int(w) for w in _DATA["two_layer"].keys())


def _tbl(wl: int):
    return _DATA["two_layer"][str(int(wl))]


def _nearest_wl(wl: int) -> int:
    return min(WAVELENGTHS, key=lambda w: abs(w - int(wl)))


def _interp(sds_mm: float, table: dict, key: str) -> float:
    # keys may be "25", "25.0", or "38"; look each up robustly
    xs = sorted(float(s) for s in table.keys())
    ys = []
    for x in xs:
        d = table.get(f"{x}") or table.get(f"{x:g}") or table.get(str(int(x)))
        ys.append(float(d[key]))
    xc = float(np.clip(sds_mm, xs[0], xs[-1]))
    return float(np.interp(xc, xs, ys))


def f_cortex_2L(sds_mm: float, wavelength_nm: int) -> float:
    """Two-layer cortical sensitivity fraction f_cortex(SDS, lambda) (production MC)."""
    wl = _nearest_wl(wavelength_nm)
    return _interp(sds_mm, _tbl(wl), "f_cortex")


def gamma_csf(sds_mm: float, wavelength_nm: int) -> float:
    """Wavelength-specific CSF light-piping ratio gamma(SDS, lambda) (production MC)."""
    wl = _nearest_wl(wavelength_nm)
    return _interp(sds_mm, _DATA["csf"][str(wl)], "gamma")


def f_cortex_invivo(sds_mm: float, wavelength_nm: int) -> float:
    """CSF-augmented (three-layer) fraction for real heads: f2L(lambda) * gamma(lambda)."""
    return min(f_cortex_2L(sds_mm, wavelength_nm) * gamma_csf(sds_mm, wavelength_nm), 0.999)


def kappa_pv(sds_mm: float, wavelength_nm: int, invivo: bool = False) -> float:
    f = f_cortex_invivo(sds_mm, wavelength_nm) if invivo else f_cortex_2L(sds_mm, wavelength_nm)
    return 1.0 / f if f > 0 else float("inf")


def provenance() -> dict:
    m = _DATA["_meta"]
    return {k: m[k] for k in ("version", "N_per_config", "n_batches", "L_max",
                              "z_max", "half_width", "g") if k in m}
