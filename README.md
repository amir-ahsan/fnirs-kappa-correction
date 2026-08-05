# Reproducibility Package — fNIRS Partial-Volume (κ) Correction

**Partial-volume correction in continuous-wave fNIRS: operating regime, a
Monte-Carlo convergence caveat, and a reproducible pipeline**

Neth Sagara and Amir Ahsan
Department of Physics, Irvine Valley College, Irvine, CA 92618, USA

This package accompanies the arXiv preprint. It contains everything needed to
(a) rebuild the manuscript PDF from source and (b) regenerate every quantitative
result, table, and figure in the paper.

---

## Overview

The study develops a transparent multiplicative correction for the systematic
bias in Modified Beer–Lambert Law (MBLL) fNIRS processing:

```
applied correction = kappa_DPF x kappa_PV        (kappa_SSR is a DIAGNOSTIC, not applied)
```

- **kappa_PV = 1 / f_cortex** corrects partial-volume dilution and is the dominant,
  and in practice the only substantial, applied factor (~9–23x).
- **kappa_DPF** corrects differential-pathlength-factor mismatch (= 1 by construction
  in the synthetic validation, which uses the true DPF in both the forward model
  and the inversion).
- **kappa_SSR = 1 / (1 − R²_SS)** is computed and **reported** per wavelength as a
  diagnostic of residual superficial contamination, but is **not applied**, because
  applying it on top of the (correct, large) kappa_PV would double-count the dilution.

Two methodological cautions are the core contribution:

1. **Numerical:** `f_cortex` must come from a *converged* forward model. The analytical
   two-layer Kienle Jacobian evaluated by fixed-point Hankel transform is **not**
   quadrature-converged for this integral and overstates `f_cortex` by ~10x. The
   package therefore calibrates `f_cortex` with a converged Monte Carlo
   (`code/mc_2layer.py`), cross-checked against a finite-difference solver and the
   published Colin27 atlas Monte Carlo.
2. **Statistical:** `kappa_SSR` over-restores and belongs in a report, not in the
   correction.

A practical consequence: the correction is reliable **only at long source–detector
separations (≥ 38 mm)**. At short separations the diluted cortical signal falls
below the measurement-noise floor and kappa_PV amplifies noise.

---

## Package structure

```
.
├── README.md                     This file
├── REVIEW_NOTES.md               Pre-submission review summary + the one edit applied
├── requirements.txt              Python dependencies
│
├── arxiv_submission/             <-- upload the CONTENTS of this folder to arXiv
│   ├── main.tex                  Manuscript source (standard article class, 11pt)
│   └── figures/
│       ├── figure1_timeseries.png
│       ├── figure2_summary.png
│       ├── figure3_robustness.png
│       ├── figure4_realdata_timeseries.png
│       ├── figure5_realdata_hrf.png
│       ├── figure6_realdata_summary.png
│       └── figure7_group_block_average.png
│
├── code/
│   ├── fnirs_kappa_synthetic_validation.py   Synthetic validation (Tables of §Results;
│   │                                         Figs 1–3; SNR, robustness, convergence,
│   │                                         CSF, model-mismatch analyses)
│   ├── fnirs_kappa_realdata_analysis.py      Single-subject real-data pipeline
│   │                                         (Subject 01; Figs 4–6; Table 6)
│   ├── fnirs_kappa_group_analysis.py         Five-subject group analysis
│   │                                         (Fig 7; Table 7; group_summary.{csv,json})
│   ├── fnirs_invivo_demo.py                  Single-channel pipeline walkthrough
│   ├── mc_2layer.py                          Converged two-layer MC → f_cortex, kappa_PV
│   └── mc_csf.py                             White MC for the CSF light-piping ratio γ
│
└── supplementary/
    ├── fnirs_kappa_beginner_notebook.ipynb  Annotated Jupyter walkthrough
    └── fNIRS_Kappa_Pedagogical_Guide.tex    LaTeX pedagogical companion
```

---

## Building the manuscript (arXiv)

The `arxiv_submission/` folder is self-contained: `main.tex` uses only standard
packages (`amsmath`, `graphicx`, `booktabs`, `multirow`, `array`, `enumitem`,
`microtype`, `hyperref`, `authblk`, `geometry`) and a bundled `figures/` directory.
The bibliography is inline (`thebibliography`), so **no `.bib`/BibTeX step is needed**.

```bash
cd arxiv_submission
pdflatex main.tex
pdflatex main.tex        # second pass resolves cross-references
```

This produces `main.pdf` (43 pages). To submit, upload the contents of
`arxiv_submission/` (i.e. `main.tex` and the `figures/` folder) as the arXiv
source. Suggested categories: **physics.med-ph** (primary), cross-list
**physics.optics** and **q-bio.NC**.

---

## Reproducing the results

### Prerequisites

- Python ≥ 3.9 (developed on 3.10; NumPy ≥ 1.22, works with NumPy 2.x)
- A LaTeX distribution (TeX Live) to build the PDF

### Install Python dependencies

```bash
pip install -r requirements.txt
```

### 1. Synthetic validation — the core quantitative results

```bash
cd code
python fnirs_kappa_synthetic_validation.py
```

Regenerates the synthetic-validation tables and Figures 1–3. The console output
reproduces (verified):

| SDS (mm) | f_cortex | kappa_PV | κ_SSR(760) | κ_SSR(850) | RMSE_MBLL | RMSE_corr |
|---------:|---------:|---------:|-----------:|-----------:|----------:|----------:|
| 25 | 0.043 | 23.42 | 1.17 |  9.82 | 2.105 | 3.270 |
| 30 | 0.064 | 15.54 | 1.19 | 10.15 | 2.059 | 2.002 |
| 35 | 0.091 | 11.02 | 1.22 | 10.52 | 1.992 | 1.228 |
| 38 | 0.102 |  9.76 | 1.22 | 10.64 | 1.977 | 1.102 |
| 40 | 0.110 |  9.12 | 1.21 | 10.52 | 1.966 | 1.087 |

Overall HbO₂ RMSE 2.020 → 1.929 µM; kappa_PV = 13.77 ± 5.43;
κ_SSR(760) = 1.199 ± 0.054, κ_SSR(850) = 10.330 ± 1.568.
HbR overall RMSE 0.653 → 0.183 µM (72.0%). Per-subject MAE 1.71 ± 0.35 µM.

**Runtime:** the core pipeline (main tables) finishes in well under a minute; the
full script also runs the grid-convergence, finite-difference, robustness, and
figure-generation steps (which use the slower analytical kernel), so end-to-end it
takes roughly **15–30 minutes** on a laptop. Set `MPLBACKEND=Agg` to run headless.

### 2. Single-subject real-data analysis (Subject 01)

```bash
cd code
python fnirs_kappa_realdata_analysis.py
```

Downloads the BIDS-NIRS-Tapping dataset (~50 MB, cached; or place a local copy at
`code/BIDS-NIRS-Tapping-data/` containing `sub-01 … sub-05` to run offline), then
applies OD conversion, 0.01–0.1 Hz bandpass, per-wavelength SSR, the kappa_SSR
diagnostic, kappa_PV correction, and MBLL inversion. Generates Figures 4–6.

Expected (SDS ≈ 38 mm): kappa_PV = 6.49 (CSF-augmented three-layer f_cortex ≈ 0.15);
κ_SSR(760) = 2.433, κ_SSR(850) = 3.303 (diagnostic only); HbO₂ 0.265 → 0.946 µM
(3.58×); HbR −0.058 → −0.181 µM (3.15×). **Runtime:** ~2–5 min (incl. first download).

Peak amplitudes are extracted within a physiological response window
(`PEAK_WINDOW_S = 2–15 s` after onset) via `windowed_peak_abs()`, rather than over the
whole epoch; for this dataset the values are identical to an unwindowed maximum.

### 3. Upgraded in-vivo analysis (primary real-data result)

```bash
cd code
python fnirs_kappa_realdata_v2.py     # per-channel, contralateral, QC'd, window-mean
```

The upgraded pipeline uses per-channel SDS/f_cortex, TDDR + SCI quality control,
the nearest short-channel regressor, condition-resolved **contralateral** channel
selection, and a fixed pre-registered window-mean estimator. It regenerates
Figures 4–7 and writes `realdata_v2_summary.json`. Expected group means (N = 5,
SDS ≈ 38 mm): corrected contralateral HbO₂ **0.84 ± 0.31 µM** (~7×), HbR
−0.05 → −0.32 µM, per-channel kappa_PV ≈ 6.3, kappa_SSR(760/850) = 1.78/2.34
(reported as a variance diagnostic, not applied). **Runtime:** ~2–4 min.

The earlier median-SDS scripts (`fnirs_kappa_realdata_analysis.py`,
`fnirs_kappa_group_analysis.py`) are retained for reference.

### 4. Monte-Carlo uncertainty and CSF ratio

```bash
cd code
python mc_uncertainty.py -N 150000 --seeds 3 --csf 1 2 3   # multi-seed f_cortex + gamma
```

Multi-seed anisotropic (g = 0.9) Monte Carlo for both wavelengths and both
geometries; reports f_cortex mean ± SD, 95% CI, effective sample size N_eff, and
per-SDS DPF, plus the CSF ratio γ at several thicknesses. Writes raw JSON/CSV
(see `results/`). **Runtime:** ~10–25 min depending on N/seeds.

### 5. Single-channel demonstration

```bash
cd code
python fnirs_invivo_demo.py
```

NumPy/Matplotlib-only walkthrough of the pipeline stages at 38 mm (kappa_PV ≈ 9.8).
**Runtime:** a few seconds.

### Optional: recompute the Monte-Carlo `f_cortex` and CSF ratio

```bash
cd code
python mc_2layer.py                 # two-layer f_cortex, kappa_PV per SDS (defaults: N=3e5, seed=1, g=0.9)
python mc_2layer.py --N 900000      # paper photon count (lower MC noise)
python mc_csf.py --geo 3L           # CSF (three-layer) f_cortex; ratio γ = f3L/f2L is the light-piping factor
python mc_2layer.py --help          # full options
```

The converged `f_cortex` values these produce are embedded as constants in the
analysis scripts (`_F_CORTEX_MC`, `_CSF_GAMMA_MC`); rerunning the MC regenerates
them (seed-insensitive; larger `--N` reduces MC noise). **Runtime:** seconds to
a couple of minutes depending on `--N`.

---

## Random seeds

- Synthetic data: `seed = 42`.
- Monte-Carlo forward models: `seed = 1` (N = 9×10⁵ photons/SDS for the two-layer
  model, 2.2×10⁶ for the CSF model in the paper).

All reported synthetic numbers are exactly reproducible with these seeds.

---

## Data and licensing

- Experimental data: BIDS-NIRS-Tapping (Luke et al., 2021),
  DOI:10.5281/zenodo.5529797 — publicly available; the real-data script fetches it
  automatically.
- Add your preferred code/data license here before public deposit (e.g. MIT for
  code, CC-BY for text/figures). No license file is included in this package.
