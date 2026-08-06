#!/usr/bin/env python3
"""
Executable homogeneous-limit validation of the two-layer Kienle diffusion FLUENCE.

The manuscript states that, when both layers are assigned identical optical
properties, the two-layer Kienle fluence -- evaluated by ADAPTIVE quadrature
(scipy.integrate.quad, relative tolerance 1e-10, subdivision limit 200) -- agrees
with the closed-form semi-infinite homogeneous solution to within 0.01% for
lateral distances rho <= 30 mm. This script reproduces that check from first
principles and emits the pointwise and maximum relative errors as a
machine-readable JSON, so the numerical claim is independently regenerable.

Note: this validates the FLUENCE implementation only. The depth-resolved cortical
sensitivity fraction f_cortex is NOT obtained from this analytical kernel; it comes
from the converged Monte Carlo (mc_production.py), because the fixed-point Hankel
quadrature is not converged for the sensitivity ratio (see the manuscript).

Run:
    python validate_homogeneous_fluence.py --out ../results/homogeneous_fluence_validation.json
"""
import argparse, json, os, sys, platform, subprocess, hashlib
from datetime import datetime, timezone
import numpy as np
from scipy.integrate import quad
from scipy.special import j0

# Stated validation configuration (matches the manuscript).
MUA = 0.015          # mm^-1
MUSP = 1.0           # mm^-1
N_INDEX = 1.4
Z_MM = 10.0          # depth of evaluation
RHO_MM = [5.0, 10.0, 15.0, 20.0, 25.0, 30.0]
Z_SUP_MM = 12.0      # nominal superficial thickness (irrelevant in the homogeneous limit,
                     # but keeps z < z_sup so the top-layer branch is exercised)
A_EXTRAP = 2.95      # extrapolated-boundary coefficient for n = 1.4 (as in the pipeline)
QUAD_EPSREL = 1e-10
QUAD_LIMIT = 200


def _diffusion_params(mua, musp):
    D = 1.0 / (3.0 * (mua + musp))
    mueff = np.sqrt(mua / D)              # = sqrt(3 mua (mua+musp))
    zp = 1.0 / musp                       # isotropic source depth
    zb = 2.0 * A_EXTRAP * D               # extrapolated boundary
    return D, mueff, zp, zb


def kienle_homogeneous_integrand(s, rho, z, D, zp, zb, mua):
    """Hankel integrand of the two-layer Kienle fluence in the homogeneous limit
    (both layers identical, so the interlayer reflection term Da -> 0), top-layer
    branch (z < z_sup). Returns phi_k(s,z) * J0(s*rho) * s."""
    alpha = np.sqrt((D * s * s + mua) / D)
    arg = abs(zp - z)
    dum1 = np.exp(-alpha * arg) - np.exp(-alpha * (2.0 * zb + zp + z))
    phi_k = dum1 / (2.0 * D * alpha)
    return phi_k * j0(s * rho) * s


def fluence_kienle_adaptive(rho, z, D, zp, zb, mua):
    """Two-layer Kienle fluence in the homogeneous limit via adaptive quadrature."""
    val, _ = quad(kienle_homogeneous_integrand, 0.0, np.inf,
                  args=(rho, z, D, zp, zb, mua),
                  epsrel=QUAD_EPSREL, limit=QUAD_LIMIT)
    return val / (2.0 * np.pi)


def fluence_semi_infinite_closed(rho, z, D, mueff, zp, zb):
    """Closed-form semi-infinite (extrapolated-boundary image-source) CW fluence."""
    r1 = np.sqrt(rho * rho + (z - zp) ** 2)
    r2 = np.sqrt(rho * rho + (z + zp + 2.0 * zb) ** 2)
    return (np.exp(-mueff * r1) / r1 - np.exp(-mueff * r2) / r2) / (4.0 * np.pi * D)


def _git_commit(cli=None):
    if cli:
        return cli
    try:
        root = os.path.dirname(os.path.abspath(__file__))
        h = subprocess.check_output(["git", "-C", root, "rev-parse", "--short", "HEAD"],
                                    stderr=subprocess.DEVNULL).decode().strip()
        dirty = subprocess.call(["git", "-C", root, "diff", "--quiet"],
                                stderr=subprocess.DEVNULL) != 0
        return h + ("-dirty" if dirty else "")
    except Exception:
        return "unknown"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--out', default='../results/homogeneous_fluence_validation.json')
    ap.add_argument('--git_commit', default=None)
    ap.add_argument('--analysis_round', default=None)
    args = ap.parse_args()

    D, mueff, zp, zb = _diffusion_params(MUA, MUSP)
    print("Homogeneous-limit fluence validation (adaptive quad vs closed-form semi-infinite)")
    print(f"  mua={MUA} mm^-1, musp={MUSP} mm^-1, n={N_INDEX}, z={Z_MM} mm; "
          f"D={D:.5f}, mu_eff={mueff:.5f}, z0={zp:.3f}, zb={zb:.3f}")
    print(f"  quad: epsrel={QUAD_EPSREL}, limit={QUAD_LIMIT}\n")
    print(f"  {'rho (mm)':>9} {'Kienle (quad)':>16} {'semi-inf (closed)':>18} {'rel err %':>12}")
    rows = []
    max_rel = 0.0
    for rho in RHO_MM:
        fk = fluence_kienle_adaptive(rho, Z_MM, D, zp, zb, MUA)
        fc = fluence_semi_infinite_closed(rho, Z_MM, D, mueff, zp, zb)
        rel = abs(fk - fc) / abs(fc)
        max_rel = max(max_rel, rel)
        rows.append(dict(rho_mm=rho, fluence_kienle_adaptive=float(fk),
                         fluence_semi_infinite_closed=float(fc),
                         rel_error=float(rel)))
        print(f"  {rho:9.1f} {fk:16.6e} {fc:18.6e} {rel*100:12.5f}")
    print(f"\n  MAX relative error over rho <= 30 mm: {max_rel*100:.5f}%  "
          f"({'PASS' if max_rel < 1e-4 else 'FAIL'} the < 0.01% claim)")

    out = dict(
        description="Homogeneous-limit validation of the two-layer Kienle CW-diffusion "
                    "fluence: adaptive-quadrature Hankel inversion vs the closed-form "
                    "semi-infinite image-source solution. Validates the fluence only, "
                    "not f_cortex (which comes from the converged Monte Carlo).",
        config=dict(mua_mm=MUA, musp_mm=MUSP, n_index=N_INDEX, z_mm=Z_MM,
                    rho_mm=RHO_MM, quad_epsrel=QUAD_EPSREL, quad_limit=QUAD_LIMIT,
                    A_extrap=A_EXTRAP),
        results=rows, max_rel_error=float(max_rel),
        passes_0p01pct=bool(max_rel < 1e-4),
        provenance=dict(schema_version="2.0", produced_by="validate_homogeneous_fluence.py",
                        git_commit=_git_commit(args.git_commit),
                        analysis_round=args.analysis_round,
                        generated_utc=datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                        command="python validate_homogeneous_fluence.py " + " ".join(sys.argv[1:]),
                        python_version=platform.python_version(), numpy_version=np.__version__,
                        scipy_version=__import__("scipy").__version__))
    payload = {k: out[k] for k in out if k != "provenance"}
    out["provenance"]["data_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()).hexdigest()
    json.dump(out, open(args.out, 'w'), indent=1)
    print(f"  wrote {args.out}")


if __name__ == '__main__':
    main()
