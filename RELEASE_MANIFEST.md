# Release Manifest

Production data provenance (from `results/fcortex_production.json`):

- `schema_version`: 2.0
- `git_commit` (data): round5
- `generated_utc`: 2026-08-06T07:31:50Z
- `data_sha256` (numeric payload): `3c34cbc66dfdf91e332c6d52dfccd8886f53f56929348b75aa56b27b4ca66217`
- `N_per_config`: 2000000, `n_batches`: 16, `seed`: 1, `g`: 0.9, `L_max`: 1200.0, `z_max`: 150.0
- batching: launch-defined (batch id = launch index mod n_batches), so batch b is the same launched-photon cohort across geometries (paired CSF-ratio / z_max stats)
- `python_version`: 3.11.15, `numpy_version`: 2.4.4

## SHA-256 of tracked files

| SHA-256 | Size (bytes) | File |
|---------|-------------:|------|
| `0bcc9c77d89fcdbf3b849ad8832b283a4519c4c5c6bfaf08c48c97723d3abc5f` | 484 | `.gitignore` |
| `293a1cd07a5479736b446494793ab37545a2028048f353835b9b4a3d88ecc55f` | 595 | `CITATION.cff` |
| `9d020d65a2756c2fc42cbbfdd7d88151cefcc1407cd1dad6814e73f5807f7d05` | 1083 | `LICENSE` |
| `0050fc32fc270930da596b499b58e5bda4da0f3a2ce2827e95938a7e52d99672` | 953 | `LICENSE-manuscript.md` |
| `2f99df1afec98b11bc9c25e9070db57ff21934cf5f7a2d79c86ec0792962dc06` | 14128 | `README.md` |
| `0a846e95f72295bd659335f7e7b4cc8cf4ad0f713b9ca435a794178d43213961` | 36537 | `REVISION_RESPONSE.md` |
| `c32766765ae6d6fd46fdeaf8a7663712e9c4e85f072068983b4f380d097205c7` | 6416 | `code/fcortex_source.py` |
| `19dfa80eb2eb154cff49449e7a047d7e26e73614986b7145c205c4840b260015` | 20804 | `code/fnirs_invivo_demo.py` |
| `546eb0b81f3dc9a652063a31da30f06450e0c9f58c6cddaf3b6255cb39c8b4a6` | 14625 | `code/fnirs_kappa_group_analysis.py` |
| `c274251a3df27129a515b2d6b6d2141fd3e8497206a6d7bc28d78101a6391702` | 38961 | `code/fnirs_kappa_realdata_analysis.py` |
| `6f856d1b64b678cbc4cf34a22ee4cf157d892865e6b6d6ec118f0cd82e845a07` | 18824 | `code/fnirs_kappa_realdata_v2.py` |
| `38be31a76543e8e938762a2ef38cb19216ea312b2abc1e8bc49471fe537659d7` | 101725 | `code/fnirs_kappa_synthetic_validation.py` |
| `3cbe561bdd8f46b851f2763c19a19efb17fd05eba77d9b2c82c1817275fa3b61` | 8598 | `code/mc_2layer.py` |
| `4544749aabea9c2834e72610cdeaba2195b5c11594a76ba959407b203b484dfe` | 9408 | `code/mc_csf.py` |
| `5ed0acbc60d5a3bf1648830c42e54ddac4d3592a31615c1c48188c9fa08f0425` | 20066 | `code/mc_production.py` |
| `ad65737f3414604580ce2d3fddb94730cc1b8d99e8c921e67b25bb25ded35ca8` | 9742 | `code/mc_robustness_sweeps.py` |
| `4be23d4f324bab85fbe7ec831930bfd8bfdc3926cb683d2143531eff8bd90d04` | 7603 | `code/mc_uncertainty.py` |
| `38f2ef56e46caeb5726141d7714bbc49ab2add6f744fc87be8cb888f06377990` | 568604 | `manuscript/figures/figure1_timeseries.png` |
| `20d522737804c3f4cd02741090503e53bdc9a5947f6fb9cf8a61c4a1f28a2aea` | 226267 | `manuscript/figures/figure2_summary.png` |
| `eba62cc6762def3edfffd939be86e4d363aa9298abea40a4380e921978dcc3a2` | 517182 | `manuscript/figures/figure3_robustness.png` |
| `d3dd45785c91b79c33650458f52514c1c6d7e1ee43e49ef92f901e61e0980f9a` | 65017 | `manuscript/figures/figure4_realdata_timeseries.png` |
| `6963ff896c7cce33ec98ed0af43c8dcbbb4303f7b6bda1dc077ea7c5e3b06e32` | 95896 | `manuscript/figures/figure5_realdata_hrf.png` |
| `17f8ab140ab7498b05a91df575451b40edc524f5a6d827212890d7f4503480e7` | 35090 | `manuscript/figures/figure6_realdata_summary.png` |
| `d9a7ebf763bbd0fd920a3a9cacfc734d4c0eec9d7fd8bfeec592c1d12262762d` | 179191 | `manuscript/figures/figure7_group_block_average.png` |
| `397d934235092dbf4ce777694e5568c99c656b709b49bf2e141863a3c05b2cbf` | 1958570 | `manuscript/main.pdf` |
| `7adb69658ddbf0f6b3d10005d6addd7d2c0d5e97ed7a3b1cc30f066ec2bb556b` | 169728 | `manuscript/main.tex` |
| `c59ba6fdd55619ffa278373ccfaca680045cbd78a829802a443e561c606f7f6b` | 342 | `requirements-lock.txt` |
| `90a1b0c629cceb68fb0900439a672b99a931e6f3ecba9ae2c8da9938a02eb9b0` | 1962 | `requirements.txt` |
| `09539f8306e7d9ec2726404bccd48a3f57eb355ff138b724fade2b75fdd3f15e` | 1059 | `results/fcortex_production.csv` |
| `0fb215ca20f486a0cbee23f183fbde6625bb4379d8d02a426aa91bcaf05b265c` | 18075 | `results/fcortex_production.json` |
| `53ad26f40b8ae6bad8936b2cbb5eee008e28836d3835c6347878881bb9c2bd7f` | 771 | `results/mc_uncertainty_2layer.csv` |
| `1f574475f248a4361fb7b4e1f162c8cc8d7728f40934813fb0e4fe7a9b21aca8` | 6389 | `results/mc_uncertainty_2layer_and_csf.json` |
| `a2bc50e82da26551ee4d6f8bd498be52bbb25e0bacda8fc861bf510db9673109` | 1066 | `results/multiseed_operating_regime.json` |
| `69f67647dbe504b8fe088cea1bf58549fd9c266e4dc7706209be402e58dbbf6e` | 19240 | `results/realdata_v2_summary.json` |
| `ba6a70efbaae023629272a6b638567b7958fe0f159f4df38f4917607e3779cd8` | 4579 | `results/robustness_secondary.json` |
| `b0f09d07fc4ef033d4b238f244f484a09477acb61bc6f560b8d9e3b7c95e420a` | 395493 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.pdf` |
| `97662882096e17270c9d5874dec90902c95e742a5e92ecbd6e61fffb2c8fe008` | 108427 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.tex` |
| `f679e7428c21b0fb107a36485e855fb5f499a3ba41e9de6f1d237ad4e44ab698` | 152016 | `supplementary/fnirs_kappa_beginner_notebook.ipynb` |

## Canonical regeneration commands

- Production forward model: `python code/mc_production.py -N 2000000 --batches 16 --out fcortex_production --git_commit <commit>`, then copy the JSON/CSV into `results/`.
- Secondary robustness sweeps: `python code/mc_robustness_sweeps.py -N 120000 --seed 1 --git_commit <commit> --out results/robustness_secondary.json`.

`fcortex_source.py` requires schema 2.0, re-verifies the payload SHA-256 on load, and resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`).
