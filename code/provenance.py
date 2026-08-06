#!/usr/bin/env python3
"""
Shared provenance helper used by every artifact-producing script in this package,
so all outputs carry the SAME provenance schema.

`provenance(produced_by, ...)` returns a dict with: schema version; the producing
script; the full and short Git SHA plus a clean/dirty flag; a human-readable
`analysis_round` label (from the argument or the ANALYSIS_ROUND environment
variable); a UTC timestamp; the exact command; the Python, OS and architecture;
NumPy (and, if importable, SciPy/pandas/MNE) versions; optional input-artifact
hashes; and the intended output destination. Producers add their run-specific
fields via `extra`.
"""
import os, sys, platform, subprocess
from datetime import datetime, timezone

SCHEMA_VERSION = "2.0"


def git_info():
    """Return (short_sha_with_dirty, full_sha, is_dirty) or ('unknown', None, None)."""
    try:
        root = os.path.dirname(os.path.abspath(__file__))
        full = subprocess.check_output(["git", "-C", root, "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
        short = subprocess.check_output(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                                        stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.call(["git", "-C", root, "diff", "--quiet"],
                                stderr=subprocess.DEVNULL) != 0
        return (short + ("-dirty" if dirty else ""), full, dirty)
    except Exception:
        return ("unknown", None, None)


def _pkg_versions():
    v = {}
    import numpy as _np
    v["numpy"] = _np.__version__
    for name in ("scipy", "pandas", "matplotlib", "mne", "mne_bids", "h5py"):
        try:
            v[name] = __import__(name).__version__
        except Exception:
            pass
    return v


def provenance(produced_by, analysis_round=None, git_commit=None,
               input_hashes=None, output=None, extra=None):
    short, full, dirty = git_info()
    if git_commit:                      # explicit override (e.g. a clean release SHA)
        short = git_commit
    block = dict(
        schema_version=SCHEMA_VERSION,
        produced_by=produced_by,
        git_commit=short,
        git_commit_full=full,
        git_dirty=dirty,
        analysis_round=analysis_round or os.environ.get("ANALYSIS_ROUND"),
        generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        command="python " + " ".join(sys.argv),
        python_version=platform.python_version(),
        platform=platform.platform(),
        machine=platform.machine(),
        packages=_pkg_versions(),
    )
    if input_hashes:
        block["input_hashes"] = input_hashes
    if output:
        block["output"] = output
    if extra:
        block.update(extra)
    return block


def fcortex_production_input():
    """Convenience: the hash/commit of the production forward model this run consumed."""
    try:
        import fcortex_source as fs
        p = fs.provenance()
        return dict(fcortex_production_sha256=p.get("data_sha256"),
                    fcortex_production_git=p.get("git_commit"),
                    fcortex_production_schema=p.get("schema_version"))
    except Exception:
        return {}
