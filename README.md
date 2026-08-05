# Partial-volume (κ) correction in continuous-wave fNIRS

**Operating regime, a Monte-Carlo convergence caveat, and a reproducible pipeline**

Neth Sagara and Amir Ahsan — Department of Physics, Irvine Valley College, Irvine, CA 92618, USA
Corresponding author: aahsan@ivc.edu

This repository contains the manuscript source and the complete analysis code for a
transparent, reproducible partial-volume correction for functional near-infrared
spectroscopy (fNIRS). It lets any reader rebuild the manuscript and regenerate every
quantitative result, table, and figure.

---

## What this study does

Routine Modified Beer–Lambert Law (MBLL) fNIRS processing treats the head as
homogeneous and underestimates cortical hemoglobin changes by roughly an order of
magnitude. We assemble a decomposed multiplicative correction and make two
methodological cautions concrete and reproducible:

```
applied correction = kappa_DPF x kappa_PV        (kappa_SSR is a DIAGNOSTIC, not applied)
```

- **kappa_PV = 1 / f_cortex** — corrects partial-volume dilution; the dominant, and in
  practice the only substantial, applied factor (~9–23x).
- **kappa_DPF** — corrects differential-pathlength mismatch (= 1 by construction in the
  synthetic validation, which uses the true DPF in forward model and inversion).
- **kappa_SSR = 1 / (1 − R²_SS)** — computed and **reported** per wavelength as a
  diagnostic of residual superficial contamination, but **not applied** (applying it on
  top of the correct, large kappa_PV would double-count the dilution).

Two cautions are the core contribution: (1) `f_cortex` must come from a *converged*
forward model — the analytical two-layer Kienle Jacobian by fixed-point Hankel transform
is not quadrature-converged for this integral and overstates `f_cortex` by ~10×, so we
calibrate against a converged Monte Carlo; and (2) `kappa_SSR` over-restores and belongs
in a report, not the correction. A practical consequence: the correction is reliable
**only at long source–detector separations (≥ 38 mm)**.

---

## Repository structure

```
.
├── README.md
├── LICENSE                     MIT — applies to the code (code/, supplementary notebook)
├── LICENSE-manuscript.md       CC-BY-4.0 — applies to the manuscript text and figures
├── CITATION.cff                Machine-readable citation metadata
├── requirements.txt            Python dependencies
│
├── manuscript/
│   ├── main.tex                Manuscript source (standard article class, arXiv-ready)
│   ├── main.pdf                Compiled manuscript (44 pp.)
│   └── figures/                figure1..figure7 (PNG)
│
├── code/
│   ├── fnirs_kappa_synthetic_validation.py   Synthetic validation (main tables; Figs 1–3)
│   ├── fnirs_kappa_realdata_analysis.py      Single-subject real-data pipeline (Figs 4–6)
│   ├── fnirs_kappa_group_analysis.py         Five-subject group analysis (Fig 7)
│   ├── fnirs_invivo_demo.py                  Single-channel pipeline walkthrough
│   ├── mc_2layer.py                          Converged two-layer MC → f_cortex, kappa_PV
│   └── mc_csf.py                             White MC for the CSF light-piping ratio γ
│
└── supplementary/
    ├── fnirs_kappa_beginner_notebook.ipynb   Annotated Jupyter walkthrough
    └── fNIRS_Kappa_Pedagogical_Guide.tex     LaTeX pedagogical companion
```

---

## Reproducing the results

Requires Python ≥ 3.9 (NumPy ≥ 1.22, works with NumPy 2.x) and, for the PDF, a LaTeX
distribution (TeX Live).

```bash
pip install -r requirements.txt

# 1) Synthetic validation — the core quantitative results (Tables; Figs 1–3)
cd code && python fnirs_kappa_synthetic_validation.py

# 2) Single-subject real data (Subject 01; Figs 4–6)
python fnirs_kappa_realdata_analysis.py

# 3) Five-subject group analysis (Fig 7; group_summary.{csv,json})
python fnirs_kappa_group_analysis.py
```

**Expected (verified) numbers.** Synthetic HbO₂ RMSE improvement −55%, +3%, +38%, +44%,
+45% at 25/30/35/38/40 mm; overall MBLL 2.020 → 1.929 µM; HbR 0.653 → 0.183 µM (72.0%);
kappa_PV = 13.77 ± 5.43; κ_SSR(760) = 1.199 ± 0.054, κ_SSR(850) = 10.330 ± 1.568.
Real data (SDS ≈ 38.3 mm): group-mean corrected HbO₂ **1.13 ± 0.41 µM** (4.59×), applied
kappa_PV = 6.49.

**Runtime.** The core synthetic correction pipeline runs in under a minute; the full
validation script (grid-convergence, finite-difference, robustness, figures via the
analytical kernel) takes ~15–30 min. Real-data scripts run in seconds excluding data
loading. Set `MPLBACKEND=Agg` to run headless.

**Random seeds.** Synthetic data `seed = 42`; Monte-Carlo forward models `seed = 1`
(9×10⁵ photons/SDS two-layer, 2.2×10⁶ CSF). All synthetic numbers are exactly reproducible.

### Building the manuscript

```bash
cd manuscript
pdflatex main.tex && pdflatex main.tex   # inline bibliography; no BibTeX step
```

---

## Data availability

The experimental data are the publicly available **BIDS-NIRS-Tapping** dataset
(Luke et al., 2021), DOI [10.5281/zenodo.5529797](https://doi.org/10.5281/zenodo.5529797).
`fnirs_kappa_realdata_analysis.py` downloads it automatically on first run (or place a
local copy at `code/BIDS-NIRS-Tapping-data/` containing `sub-01 … sub-05`). It is **not**
redistributed in this repository.

---

## Licensing

This repository is dual-licensed:

- **Code** (`code/`, `supplementary/*.ipynb`) — **MIT License** (see `LICENSE`).
- **Manuscript and figures** (`manuscript/`, `supplementary/*.tex`) — **Creative Commons
  Attribution 4.0 International (CC-BY-4.0)** (see `LICENSE-manuscript.md`).

You are free to reuse both with attribution.

---

## How to cite

See `CITATION.cff` (GitHub renders a "Cite this repository" button). A preprint DOI /
arXiv identifier will be added here once the preprint is posted.

## AI use disclosure

A large language model was used as an assistive tool for code development/debugging,
running and cross-checking analyses, and drafting/revising text; the authors directed all
scientific decisions and verified all code, results, and text. See the manuscript's
"Disclosure of the Use of Artificial Intelligence" section.
