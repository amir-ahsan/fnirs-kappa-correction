# Release Manifest

Production data provenance (from `results/fcortex_production.json`):

- `schema_version`: 2.0
- `git_commit` (data, real SHA): 0782d6e  |  `analysis_round`: round8
- `generated_utc`: 2026-08-07T05:16:54Z
- `data_sha256` (numeric payload): `91cfa828227df0652fb023a9ec5a005424103bcc5841c162fe9cfc58d2ac0331`
- SDS grid (mm): [25.0, 30.0, 33.0, 35.0, 38.0, 40.0, 42.0]  (spans the in-vivo channel range 33.4-40.9 mm)
- `N_per_config`: 2000000, `n_batches`: 16, `seed`: 1, `N_thin_csf`: 2000000, `g`: 0.9, `L_max`: 1200.0, `z_max`: 150.0
- `python_version`: 3.11.15  |  platform: Linux-6.18.5-fc-v18-x86_64-with-glibc2.39
- packages: numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, matplotlib 3.10.9, mne 1.12.1, mne-bids 0.19.0, h5py 3.16.0

All artifacts share the uniform provenance schema written by `code/provenance.py`.

## SHA-256 of tracked files

| SHA-256 | Size (bytes) | File |
|---------|-------------:|------|
| `0bcc9c77d89fcdbf3b849ad8832b283a4519c4c5c6bfaf08c48c97723d3abc5f` | 484 | `.gitignore` |
| `293a1cd07a5479736b446494793ab37545a2028048f353835b9b4a3d88ecc55f` | 595 | `CITATION.cff` |
| `9d020d65a2756c2fc42cbbfdd7d88151cefcc1407cd1dad6814e73f5807f7d05` | 1083 | `LICENSE` |
| `0050fc32fc270930da596b499b58e5bda4da0f3a2ce2827e95938a7e52d99672` | 953 | `LICENSE-manuscript.md` |
| `0617f8a524d2d76df5256659daf081c668add09356e2e774eefa796584911f64` | 16277 | `README.md` |
| `4e4048ed8c3ecb57208a66838f7539d9c6103720127b8524129b185f70e985b8` | 56018 | `REVISION_RESPONSE.md` |
| `1f9705d3eba85093b106f62b46e08bc7e10c9bf44c95bbc396bf85a1f39fa75b` | 6985 | `code/fcortex_source.py` |
| `19dfa80eb2eb154cff49449e7a047d7e26e73614986b7145c205c4840b260015` | 20804 | `code/fnirs_invivo_demo.py` |
| `546eb0b81f3dc9a652063a31da30f06450e0c9f58c6cddaf3b6255cb39c8b4a6` | 14625 | `code/fnirs_kappa_group_analysis.py` |
| `48c04320d47460f6a77941967d5fcdb0a8b767ed76357cd7ecf46c5b1f52af60` | 42500 | `code/fnirs_kappa_realdata_analysis.py` |
| `ca9c5572a9433b6a3d41f748b88e72af1ee513f9cc48a8c9037aa6021656e5f7` | 25660 | `code/fnirs_kappa_realdata_v2.py` |
| `f966843d5a001ba50e228ec0740d2c3365fff8838607b84d2faedf8f425469b4` | 105148 | `code/fnirs_kappa_synthetic_validation.py` |
| `3cbe561bdd8f46b851f2763c19a19efb17fd05eba77d9b2c82c1817275fa3b61` | 8598 | `code/mc_2layer.py` |
| `4544749aabea9c2834e72610cdeaba2195b5c11594a76ba959407b203b484dfe` | 9408 | `code/mc_csf.py` |
| `3fd7aa0a5a2c6d73e9469153d98f4e1596eb079fb212e1eac869ee56161e08ca` | 22602 | `code/mc_production.py` |
| `55ec336ad7ff21e23194e445dbd7ab55a2b33c19f8e40d8db5844dc7a0a7c09c` | 9843 | `code/mc_robustness_sweeps.py` |
| `4be23d4f324bab85fbe7ec831930bfd8bfdc3926cb683d2143531eff8bd90d04` | 7603 | `code/mc_uncertainty.py` |
| `4488f12aaa353f2e2787cf1e15fe31a73e909a83ab23775c477cc3425dc6f475` | 3252 | `code/provenance.py` |
| `f5c75ecf36755d119c93ddf523555c372c4cce46548013c219935a3f6d98ae7c` | 8290 | `code/validate_homogeneous_fluence.py` |
| `38f2ef56e46caeb5726141d7714bbc49ab2add6f744fc87be8cb888f06377990` | 568604 | `manuscript/figures/figure1_timeseries.png` |
| `20d522737804c3f4cd02741090503e53bdc9a5947f6fb9cf8a61c4a1f28a2aea` | 226267 | `manuscript/figures/figure2_summary.png` |
| `95b6b39f7350a1bd8eedcb235ddc80d51f703116f391d75614bec13a57689ad6` | 533977 | `manuscript/figures/figure3_robustness.png` |
| `f98aea76b972f92984f2e7f27f39962149938bc1b8c5219310164586b17c8bbb` | 64503 | `manuscript/figures/figure4_realdata_timeseries.png` |
| `969c0176606c8700ff3353f24a83c9ac59b28a5f8169ef2492875067e46f2af8` | 96194 | `manuscript/figures/figure5_realdata_hrf.png` |
| `8eee074082111c72b5f38b02b3668209c37a488afbdebe51276387c7c6fcfc22` | 35028 | `manuscript/figures/figure6_realdata_summary.png` |
| `d0207285253801dd05350a8bb54700492a5428bf3b3c2ff8557278e13791954a` | 178510 | `manuscript/figures/figure7_group_block_average.png` |
| `c661b477d5bb3833b35d4802d0a0dacc184cd57e971ece582b45208ae120f4ac` | 1980282 | `manuscript/main.pdf` |
| `432d522d109ab4acf916b6e9c9c0eef33ed519550b948f821d57c3c09af8a439` | 176354 | `manuscript/main.tex` |
| `df231cbbdc86ea041702c58bc98ca791be7a2d12640fbb298308ce600a349dbb` | 4575 | `reproduce_all.sh` |
| `f50e014ad9afaa7cd8a0d1c79b91d516570ac87e72ce4c8d052f862ba7d793a1` | 2687 | `requirements-lock.txt` |
| `90a1b0c629cceb68fb0900439a672b99a931e6f3ecba9ae2c8da9938a02eb9b0` | 1962 | `requirements.txt` |
| `2e06607221e8e81db393b5eb73f915f6c69e9ef95ba2844ab8f60441472fa0db` | 2701 | `results/fcortex_production.csv` |
| `d8e87af5e584b63d98c71e812bcbbad714cca48572951964394d1a0a9f983e09` | 30630 | `results/fcortex_production.json` |
| `fc7e2dc760cb951f84fa46b325a63fbbb09cbe29870d28019a203718e406e494` | 3687 | `results/homogeneous_fluence_validation.json` |
| `4b3424f183861894d3c966302c4baeb977e9e1097c462480c54d0cb7580e52ba` | 770 | `results/legacy/README.md` |
| `53ad26f40b8ae6bad8936b2cbb5eee008e28836d3835c6347878881bb9c2bd7f` | 771 | `results/legacy/mc_uncertainty_2layer.csv` |
| `1f574475f248a4361fb7b4e1f162c8cc8d7728f40934813fb0e4fe7a9b21aca8` | 6389 | `results/legacy/mc_uncertainty_2layer_and_csf.json` |
| `eb563320250dadf2419d7288737a4ab4f50941cea5e8b01a752b56569fd293fd` | 2154 | `results/multiseed_operating_regime.json` |
| `ac5904095c7d258183f94494246e5a93b2e534ddb86268b6405a0cc1dbbaac69` | 32743 | `results/realdata_v2_summary.json` |
| `1486c97e6caaea3ad3e8906639a3c7ebe67a36e26ed14df09724b1aa0b3cc56a` | 4987 | `results/robustness_secondary.json` |
| `b459a6309f976c07c7f9e401e9d83df92e91e9c2ae5c0bce89b8e0378ba09341` | 420754 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.pdf` |
| `81d28dc91e482cbd54822c80e75358197a73ad9051e0d0e14ff2b7cd50e04cb9` | 110112 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.tex` |
| `928e138ed7d676fca0a6ecdcd61ec4a4c3556479b0b27865a617d350b2db00df` | 152680 | `supplementary/fnirs_kappa_beginner_notebook.ipynb` |

## Reproduce everything

Run `./reproduce_all.sh` (or `FAST=1 ./reproduce_all.sh` for a smoke run). Canonical commands:

- Forward model: `python code/mc_production.py -N 2000000 --batches 16 --thin_csf_N 2000000 --out fcortex_production --git_commit <sha> --analysis_round <label>` → copy JSON/CSV to `results/`.
- Robustness (Fig 3 source): `python code/mc_robustness_sweeps.py -N 500000 --seed 1 --git_commit <sha> --analysis_round <label> --out results/robustness_secondary.json`.
- Fluence validation: `python code/validate_homogeneous_fluence.py --out results/homogeneous_fluence_validation.json`.

`fcortex_source.py` requires schema 2.0, requires and re-verifies the payload SHA-256, and resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`).
