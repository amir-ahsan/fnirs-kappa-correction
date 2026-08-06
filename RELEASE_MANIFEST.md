# Release Manifest

Production data provenance (from `results/fcortex_production.json`):

- `schema_version`: 2.0
- `git_commit` (data, real SHA): 4ea6693
- `analysis_round`: round6
- `generated_utc`: 2026-08-06T18:25:12Z
- `data_sha256` (numeric payload): `e2e87a35e10f7252ecf7f10e6e6f279116ce67a367f57210db9b3ab7fc0d7bd2`
- `N_per_config`: 2000000, `n_batches`: 16, `seed`: 1, `N_thin_csf`: 2000000, `g`: 0.9, `L_max`: 1200.0, `z_max`: 150.0
- batching: launch-defined (batch id = launch index mod n_batches); the independent-propagation γ SE is primary, the launch-index-matched batch SE is a cross-check
- `python_version`: 3.11.15, `numpy_version`: 2.4.4

## SHA-256 of tracked files

| SHA-256 | Size (bytes) | File |
|---------|-------------:|------|
| `0bcc9c77d89fcdbf3b849ad8832b283a4519c4c5c6bfaf08c48c97723d3abc5f` | 484 | `.gitignore` |
| `293a1cd07a5479736b446494793ab37545a2028048f353835b9b4a3d88ecc55f` | 595 | `CITATION.cff` |
| `9d020d65a2756c2fc42cbbfdd7d88151cefcc1407cd1dad6814e73f5807f7d05` | 1083 | `LICENSE` |
| `0050fc32fc270930da596b499b58e5bda4da0f3a2ce2827e95938a7e52d99672` | 953 | `LICENSE-manuscript.md` |
| `eb9b8b579ce34f35a6f8428eb220c2b2c2b50935f369b16e0df3a4c78582651e` | 16274 | `README.md` |
| `d5ca91607c81fbda0f9b7b3e1226f2c4952da41602c3e9f69a902afd4097573c` | 43423 | `REVISION_RESPONSE.md` |
| `1f9705d3eba85093b106f62b46e08bc7e10c9bf44c95bbc396bf85a1f39fa75b` | 6985 | `code/fcortex_source.py` |
| `19dfa80eb2eb154cff49449e7a047d7e26e73614986b7145c205c4840b260015` | 20804 | `code/fnirs_invivo_demo.py` |
| `546eb0b81f3dc9a652063a31da30f06450e0c9f58c6cddaf3b6255cb39c8b4a6` | 14625 | `code/fnirs_kappa_group_analysis.py` |
| `48c04320d47460f6a77941967d5fcdb0a8b767ed76357cd7ecf46c5b1f52af60` | 42500 | `code/fnirs_kappa_realdata_analysis.py` |
| `38612f9b9a526c1173b6a4492ab3a397007f44ea8cc1d9138fcf9cc09652d509` | 20861 | `code/fnirs_kappa_realdata_v2.py` |
| `fe89407c1d00de3941fd77b2e747194ab8d1377d291801bedafa26916798ec22` | 106287 | `code/fnirs_kappa_synthetic_validation.py` |
| `3cbe561bdd8f46b851f2763c19a19efb17fd05eba77d9b2c82c1817275fa3b61` | 8598 | `code/mc_2layer.py` |
| `4544749aabea9c2834e72610cdeaba2195b5c11594a76ba959407b203b484dfe` | 9408 | `code/mc_csf.py` |
| `4de48d89ef1453daa72953b70a644660bfcc9374abff736ff7f6f4f1cdbe49ec` | 22401 | `code/mc_production.py` |
| `31090cf21da27124343e2875cc02b1779da8b177a53eea76c9ad2df395420a85` | 10017 | `code/mc_robustness_sweeps.py` |
| `4be23d4f324bab85fbe7ec831930bfd8bfdc3926cb683d2143531eff8bd90d04` | 7603 | `code/mc_uncertainty.py` |
| `701396f88c0a94141c086d95b53fd54226bc580c6aa07f5b30be1a3123b09bff` | 6835 | `code/validate_homogeneous_fluence.py` |
| `38f2ef56e46caeb5726141d7714bbc49ab2add6f744fc87be8cb888f06377990` | 568604 | `manuscript/figures/figure1_timeseries.png` |
| `20d522737804c3f4cd02741090503e53bdc9a5947f6fb9cf8a61c4a1f28a2aea` | 226267 | `manuscript/figures/figure2_summary.png` |
| `95b6b39f7350a1bd8eedcb235ddc80d51f703116f391d75614bec13a57689ad6` | 533977 | `manuscript/figures/figure3_robustness.png` |
| `d3dd45785c91b79c33650458f52514c1c6d7e1ee43e49ef92f901e61e0980f9a` | 65017 | `manuscript/figures/figure4_realdata_timeseries.png` |
| `6963ff896c7cce33ec98ed0af43c8dcbbb4303f7b6bda1dc077ea7c5e3b06e32` | 95896 | `manuscript/figures/figure5_realdata_hrf.png` |
| `17f8ab140ab7498b05a91df575451b40edc524f5a6d827212890d7f4503480e7` | 35090 | `manuscript/figures/figure6_realdata_summary.png` |
| `d9a7ebf763bbd0fd920a3a9cacfc734d4c0eec9d7fd8bfeec592c1d12262762d` | 179191 | `manuscript/figures/figure7_group_block_average.png` |
| `723bda1c59854383fc94045e80a7b5029b22d4d8eece8513fe264e0366e5041c` | 1975957 | `manuscript/main.pdf` |
| `325acbbce17064513ec783852b5729d333a681b682031cafb5509f3856bfc99a` | 172412 | `manuscript/main.tex` |
| `c59ba6fdd55619ffa278373ccfaca680045cbd78a829802a443e561c606f7f6b` | 342 | `requirements-lock.txt` |
| `90a1b0c629cceb68fb0900439a672b99a931e6f3ecba9ae2c8da9938a02eb9b0` | 1962 | `requirements.txt` |
| `f376cb22af11b572ed82ac0567ab636bfc5ccb97fd1a8e6b801e61b6586f0bc7` | 1969 | `results/fcortex_production.csv` |
| `06cd27f00dffee9e9b09b7af8be12d372b216631a1caeb31b5cc5a126d3db147` | 22335 | `results/fcortex_production.json` |
| `1f61e47d2acb92f1c4310d1b27a69155310398f7b6f712c79ebca5e8d4d80a74` | 2170 | `results/homogeneous_fluence_validation.json` |
| `4b3424f183861894d3c966302c4baeb977e9e1097c462480c54d0cb7580e52ba` | 770 | `results/legacy/README.md` |
| `53ad26f40b8ae6bad8936b2cbb5eee008e28836d3835c6347878881bb9c2bd7f` | 771 | `results/legacy/mc_uncertainty_2layer.csv` |
| `1f574475f248a4361fb7b4e1f162c8cc8d7728f40934813fb0e4fe7a9b21aca8` | 6389 | `results/legacy/mc_uncertainty_2layer_and_csf.json` |
| `5b3b5424c4981c8aadd7ab0e5f93aeb673dfd7a73ce9914c7c24ff95d4fba3a8` | 1633 | `results/multiseed_operating_regime.json` |
| `b9f80d019824984123fbf8decb2652ec68743afa96400476d13287e34ea83569` | 20268 | `results/realdata_v2_summary.json` |
| `4b4ec1eb9702a1359a4df155e4769ed4e5d37f4b34812baf3b0cd61cae9b0a2f` | 4639 | `results/robustness_secondary.json` |
| `84c1623835d03efc48b21e532a63ea908a75e2ba5e48317f54ca697f1a9f6865` | 395807 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.pdf` |
| `2bbca4aa2b28f75a09ba812aa7183b533ccc494ba20f860ea9b85e488f5cf87a` | 108761 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.tex` |
| `f679e7428c21b0fb107a36485e855fb5f499a3ba41e9de6f1d237ad4e44ab698` | 152016 | `supplementary/fnirs_kappa_beginner_notebook.ipynb` |

## Canonical regeneration commands

- Forward model: `python code/mc_production.py -N 2000000 --batches 16 --thin_csf_N 2000000 --out fcortex_production --git_commit <sha> --analysis_round <label>`, then copy JSON/CSV into `results/`.
- Robustness sweeps (Figure 3 source): `python code/mc_robustness_sweeps.py -N 500000 --seed 1 --git_commit <sha> --analysis_round <label> --out results/robustness_secondary.json`.
- Fluence validation: `python code/validate_homogeneous_fluence.py --out results/homogeneous_fluence_validation.json`.

`fcortex_source.py` requires schema 2.0, requires and re-verifies the payload SHA-256 on load, and resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`).
