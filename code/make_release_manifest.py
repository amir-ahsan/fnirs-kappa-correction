#!/usr/bin/env python3
"""
Regenerate RELEASE_MANIFEST.md from the current tracked files.

The manifest header echoes the production forward-model provenance
(`results/fcortex_production.json`) and the body lists the SHA-256 and size of
every Git-tracked file (except the manifest itself). Running this after a release
build guarantees the manifest reflects exactly what is committed, so a stale
manifest cannot survive a release. Invoked automatically by `reproduce_all.sh`
in RELEASE mode, or run directly:

    python code/make_release_manifest.py            # writes ./RELEASE_MANIFEST.md
"""
import json, hashlib, subprocess, os, sys


def repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        return subprocess.check_output(["git", "-C", here, "rev-parse", "--show-toplevel"],
                                       stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return os.path.dirname(here)


def main():
    root = repo_root()
    prod = os.path.join(root, "results", "fcortex_production.json")
    d = json.load(open(prod))["_meta"]
    pk = d.get("packages") or {}
    files = subprocess.check_output(["git", "-C", root, "ls-files"]).decode().split()
    files = sorted(f for f in files if f != "RELEASE_MANIFEST.md")

    L = ["# Release Manifest\n",
         "Production data provenance (from `results/fcortex_production.json`):\n",
         f"- `schema_version`: {d.get('data_version', '2.0')}",
         f"- `git_commit` (data, real SHA): {d.get('git_commit')}  |  `git_commit_full`: "
         f"{d.get('git_commit_full')}  |  `git_dirty`: {d.get('git_dirty')}",
         f"- `analysis_round`: {d.get('analysis_round')}"
         + (f"  |  `release_label`: {d.get('release_label')}" if d.get('release_label') else ""),
         f"- `generated_utc`: {d.get('generated_utc')}",
         f"- `data_sha256` (numeric payload): `{d.get('data_sha256')}`",
         f"- SDS grid (mm): {d.get('sds')}  (spans the in-vivo channel range 33.4-40.9 mm)",
         f"- `N_per_config`: {d.get('N_per_config')}, `n_batches`: {d.get('n_batches')}, "
         f"`seed`: {d.get('seed')}, `N_thin_csf`: {d.get('N_thin_csf')}, `g`: {d.get('g')}, "
         f"`L_max`: {d.get('L_max')}, `z_max`: {d.get('z_max')}",
         f"- `python_version`: {d.get('python_version')}  |  platform: {d.get('platform')}",
         f"- packages: numpy {pk.get('numpy')}, scipy {pk.get('scipy')}, pandas {pk.get('pandas')}, "
         f"matplotlib {pk.get('matplotlib')}, mne {pk.get('mne')}, mne-bids {pk.get('mne_bids')}, "
         f"h5py {pk.get('h5py')}",
         "",
         "All artifacts share the uniform provenance schema written by `code/provenance.py` "
         "(git_commit / git_commit_full / git_dirty are auto-detected from HEAD; a human release "
         "name is stored separately as `release_label`).\n",
         "## SHA-256 of tracked files\n",
         "| SHA-256 | Size (bytes) | File |",
         "|---------|-------------:|------|"]
    for f in files:
        fp = os.path.join(root, f)
        if not os.path.exists(fp):
            continue
        b = open(fp, "rb").read()
        L.append(f"| `{hashlib.sha256(b).hexdigest()}` | {len(b)} | `{f}` |")
    L += ["",
          "## Reproduce everything\n",
          "Run `./reproduce_all.sh` (or `FAST=1 ./reproduce_all.sh` for a smoke run; "
          "`RELEASE=1 ./reproduce_all.sh` makes the in-vivo step fatal and regenerates this "
          "manifest automatically). Canonical commands:\n",
          "- Forward model: `python code/mc_production.py -N 2000000 --batches 16 --thin_csf_N 2000000 "
          "--out fcortex_production --analysis_round <label>` → copy JSON/CSV to `results/`.",
          "- Robustness (Fig 3 source): `python code/mc_robustness_sweeps.py -N 500000 --seed 1 "
          "--analysis_round <label> --out results/robustness_secondary.json`.",
          "- Fluence validation: `python code/validate_homogeneous_fluence.py "
          "--out results/homogeneous_fluence_validation.json`.",
          "",
          "`fcortex_source.py` requires schema 2.0, requires and re-verifies the payload SHA-256, and "
          "resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`)."]
    out = os.path.join(root, "RELEASE_MANIFEST.md")
    open(out, "w").write("\n".join(L) + "\n")
    print(f"wrote {out} ({len(files)} files hashed)")


if __name__ == "__main__":
    main()
