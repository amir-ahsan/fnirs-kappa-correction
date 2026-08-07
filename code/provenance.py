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
    """Return (short_sha, full_sha, is_dirty) or ('unknown', None, None).

    The short SHA is the plain `git rev-parse --short HEAD`; the working-tree
    clean/dirty state is returned SEPARATELY as `is_dirty` (a bool), not encoded
    into the short string, so provenance records the two facts independently and
    a caller can never conflate a label with the real commit. `is_dirty` uses
    `git status --porcelain`, so it accounts for staged, unstaged AND untracked
    (non-ignored) changes -- a stricter clean-release check than a bare `git diff`
    (git-ignored files such as the downloaded dataset do not count as dirty)."""
    try:
        root = os.path.dirname(os.path.abspath(__file__))
        full = subprocess.check_output(["git", "-C", root, "rev-parse", "HEAD"],
                                       stderr=subprocess.DEVNULL).decode().strip()
        short = subprocess.check_output(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                                        stderr=subprocess.DEVNULL).decode().strip()
        # dirty = any staged, unstaged OR untracked (non-ignored) change in the tree.
        # `git status --porcelain` lists all three; an empty output means clean.
        porcelain = subprocess.check_output(["git", "-C", root, "status", "--porcelain"],
                                            stderr=subprocess.DEVNULL).decode()
        dirty = bool(porcelain.strip())
        return (short, full, dirty)
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


def provenance(produced_by, analysis_round=None, git_commit=None, release_label=None,
               input_hashes=None, output=None, extra=None):
    """Build the provenance block.

    `git_commit`, `git_commit_full` and `git_dirty` are ALWAYS determined
    automatically from `git rev-parse HEAD`; they are never overridden by a
    caller, so provenance cannot claim a commit the checkout does not actually
    have. A caller MAY pass `git_commit` as an ASSERTION that the release is being
    built from that commit: it is validated against HEAD and a mismatch RAISES
    (build from a clean checkout instead). A human-readable release name should be
    passed as `release_label` (or the RELEASE_LABEL env var) and is stored in its
    own field, distinct from the real SHA."""
    short, full, dirty = git_info()
    if git_commit and full and full != "None":
        # Treat a passed value as an assertion about HEAD; accept a full SHA, an
        # abbreviated SHA (prefix of the full), or the exact short SHA. Reject a
        # mismatch rather than silently recording a commit the tree isn't on.
        gc = str(git_commit).strip()
        if not (gc == short or full.startswith(gc) or gc.startswith(short)):
            raise ValueError(
                f"requested git_commit {gc!r} does not match HEAD ({short}/{full}); "
                "commit your changes and build the release from a clean checkout at HEAD, "
                "or pass a human label via release_label= instead of git_commit=.")
    label = release_label or os.environ.get("RELEASE_LABEL")
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
    if label:
        block["release_label"] = label
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
