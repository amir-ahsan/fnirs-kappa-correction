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
import json, os, hashlib
from pathlib import Path
import numpy as np

_PROD_NAME = "fcortex_production.json"

# This manuscript release is pinned to the schema-2.0 production artifact. Schema
# 1.0 (which lacked the per-batch paired fields and the payload hash) is rejected
# so a stale pre-2.0 file cannot silently reproduce different numbers.
REQUIRED_SCHEMA = "2.0"


def _find_production() -> Path:
    """Resolve the production file, preferring the version-controlled release
    artifact over arbitrary current-working-directory files.

    Order: an explicit path in the FCORTEX_PRODUCTION_JSON environment variable;
    then the tracked ``results/`` copy (the frozen release artifact) relative to
    this module and to the CWD; then a locally regenerated copy next to this
    module or in the CWD. Preferring ``results/`` first means a stale loose file
    in the code directory or CWD cannot shadow the released, hash-stamped data."""
    env = os.environ.get("FCORTEX_PRODUCTION_JSON")
    here = Path(__file__).resolve().parent
    candidates = []
    if env:
        candidates.append(Path(env).expanduser())
    candidates += [
        here.parent / "results" / _PROD_NAME,   # tracked release artifact (preferred)
        Path.cwd() / "results" / _PROD_NAME,
        here / _PROD_NAME,                       # locally regenerated (git-ignored)
        Path.cwd() / _PROD_NAME,
    ]
    for cand in candidates:
        if cand.exists():
            return cand
    raise FileNotFoundError(
        f"{_PROD_NAME} not found via FCORTEX_PRODUCTION_JSON, under results/, next "
        f"to fcortex_source.py, or in the CWD. Generate it first:  python "
        f"mc_production.py -N 2000000 --batches 16 "
        f"--out {Path(__file__).resolve().parent / 'fcortex_production'}")


def _payload_sha256(data: dict) -> str:
    """SHA-256 over the numeric payload only (everything except _meta), matching
    mc_production._payload_sha256, so the stored hash can be re-verified on load."""
    payload = {k: data[k] for k in data if k != "_meta"}
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(blob).hexdigest()


def _load_and_check(path: Path) -> dict:
    data = json.load(open(path))
    meta = data.get("_meta", {})
    ver = str(meta.get("schema_version", meta.get("version", "1.0")))
    if ver != REQUIRED_SCHEMA:
        raise ValueError(
            f"{path.name} has schema_version {ver!r}, but this manuscript release "
            f"requires schema {REQUIRED_SCHEMA!r}. Regenerate the production file with "
            f"the current mc_production.py (which stamps schema {REQUIRED_SCHEMA}).")
    for key in ("two_layer", "csf"):
        if key not in data:
            raise ValueError(f"{path.name} is missing the required '{key}' block "
                             f"(schema_version {ver}); it is not a valid production file.")
    # Validate the stored numeric-payload hash so an edited/corrupted/truncated
    # file cannot masquerade as the released artifact. In release mode the hash
    # is REQUIRED: a schema-2.0 file with data_sha256 removed is rejected, so the
    # integrity guard cannot be silently defeated by deleting the field. Set
    # FCORTEX_SKIP_HASH=1 for intentional hand-editing during development.
    dev_override = os.environ.get("FCORTEX_SKIP_HASH")
    stored = meta.get("data_sha256")
    if not stored and not dev_override:
        raise ValueError(
            f"{path.name} has no _meta.data_sha256 (schema {ver}). The manuscript "
            f"release requires the payload hash for integrity; regenerate the file "
            f"with mc_production.py, or set FCORTEX_SKIP_HASH=1 to override in "
            f"development.")
    if stored and not dev_override:
        actual = _payload_sha256(data)
        if actual != stored:
            raise ValueError(
                f"{path.name} payload SHA-256 mismatch: stored {stored[:12]}... but "
                f"recomputed {actual[:12]}.... The production file does not match its "
                f"provenance hash (edited, corrupted, or regenerated without updating "
                f"_meta). Regenerate it with mc_production.py or set FCORTEX_SKIP_HASH=1 "
                f"to override.")
    return data


_DATA = _load_and_check(_find_production())
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
    keys = ("schema_version", "data_version", "version", "git_commit", "analysis_round",
            "generated_utc", "command", "python_version", "numpy_version", "data_sha256",
            "N_per_config", "n_batches", "seed", "N_thin_csf", "L_max", "z_max", "half_width", "g")
    return {k: m[k] for k in keys if k in m}
