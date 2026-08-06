# Release Manifest

Production data provenance (from `results/fcortex_production.json`):

- `schema_version`: 2.0
- `git_commit` (data, real SHA): bd007a1  |  `analysis_round`: round7
- `generated_utc`: 2026-08-06T22:21:43Z
- `data_sha256` (numeric payload): `d2f2f2d6ac0e21453c2932c9109101d7af099f9325c5db1fed6dca98d1a25ad6`
- SDS grid (mm): [25.0, 30.0, 33.0, 35.0, 38.0, 40.0, 42.0]  (spans the in-vivo channel range 33.4-40.9 mm)
- `N_per_config`: 2000000, `n_batches`: 16, `seed`: 1, `N_thin_csf`: 2000000, `g`: 0.9, `L_max`: 1200.0, `z_max`: 150.0
- `python_version`: 3.11.15  |  platform: Linux-6.18.5-fc-v18-x86_64-with-glibc2.39
- packages: numpy 2.4.4, scipy 1.17.1

All artifacts share the uniform provenance schema written by `code/provenance.py`.

## SHA-256 of tracked files

| SHA-256 | Size (bytes) | File |
|---------|-------------:|------|
| `0bcc9c77d89fcdbf3b849ad8832b283a4519c4c5c6bfaf08c48c97723d3abc5f` | 484 | `.gitignore` |
| `293a1cd07a5479736b446494793ab37545a2028048f353835b9b4a3d88ecc55f` | 595 | `CITATION.cff` |
| `9d020d65a2756c2fc42cbbfdd7d88151cefcc1407cd1dad6814e73f5807f7d05` | 1083 | `LICENSE` |
| `0050fc32fc270930da596b499b58e5bda4da0f3a2ce2827e95938a7e52d99672` | 953 | `LICENSE-manuscript.md` |
| `2892f9e506b02961c516663df3c0eb59e7810a6d2d31ce967f9d0e58e8c225c7` | 16074 | `README.md` |
| `7250539b4389710aaa96b642f869739041a9c4cce08bb81f4ab1ae7c3abe0bb3` | 49989 | `REVISION_RESPONSE.md` |
| `1f9705d3eba85093b106f62b46e08bc7e10c9bf44c95bbc396bf85a1f39fa75b` | 6985 | `code/fcortex_source.py` |
| `19dfa80eb2eb154cff49449e7a047d7e26e73614986b7145c205c4840b260015` | 20804 | `code/fnirs_invivo_demo.py` |
| `546eb0b81f3dc9a652063a31da30f06450e0c9f58c6cddaf3b6255cb39c8b4a6` | 14625 | `code/fnirs_kappa_group_analysis.py` |
| `48c04320d47460f6a77941967d5fcdb0a8b767ed76357cd7ecf46c5b1f52af60` | 42500 | `code/fnirs_kappa_realdata_analysis.py` |
| `72cde10959d5271971db7c0541c0cf485e08e63e0f3617eef6f9dc3a318adcd1` | 22626 | `code/fnirs_kappa_realdata_v2.py` |
| `f966843d5a001ba50e228ec0740d2c3365fff8838607b84d2faedf8f425469b4` | 105148 | `code/fnirs_kappa_synthetic_validation.py` |
| `3cbe561bdd8f46b851f2763c19a19efb17fd05eba77d9b2c82c1817275fa3b61` | 8598 | `code/mc_2layer.py` |
| `4544749aabea9c2834e72610cdeaba2195b5c11594a76ba959407b203b484dfe` | 9408 | `code/mc_csf.py` |
| `f5c2a68b763b587481c2d2d4253cd4ec1938f5107c91797a9a4cde4d4d07a26a` | 22566 | `code/mc_production.py` |
| `55ec336ad7ff21e23194e445dbd7ab55a2b33c19f8e40d8db5844dc7a0a7c09c` | 9843 | `code/mc_robustness_sweeps.py` |
| `4be23d4f324bab85fbe7ec831930bfd8bfdc3926cb683d2143531eff8bd90d04` | 7603 | `code/mc_uncertainty.py` |
| `4488f12aaa353f2e2787cf1e15fe31a73e909a83ab23775c477cc3425dc6f475` | 3252 | `code/provenance.py` |
| `f5c75ecf36755d119c93ddf523555c372c4cce46548013c219935a3f6d98ae7c` | 8290 | `code/validate_homogeneous_fluence.py` |
| `38f2ef56e46caeb5726141d7714bbc49ab2add6f744fc87be8cb888f06377990` | 568604 | `manuscript/figures/figure1_timeseries.png` |
| `20d522737804c3f4cd02741090503e53bdc9a5947f6fb9cf8a61c4a1f28a2aea` | 226267 | `manuscript/figures/figure2_summary.png` |
| `95b6b39f7350a1bd8eedcb235ddc80d51f703116f391d75614bec13a57689ad6` | 533977 | `manuscript/figures/figure3_robustness.png` |
| `d3dd45785c91b79c33650458f52514c1c6d7e1ee43e49ef92f901e61e0980f9a` | 65017 | `manuscript/figures/figure4_realdata_timeseries.png` |
| `6963ff896c7cce33ec98ed0af43c8dcbbb4303f7b6bda1dc077ea7c5e3b06e32` | 95896 | `manuscript/figures/figure5_realdata_hrf.png` |
| `17f8ab140ab7498b05a91df575451b40edc524f5a6d827212890d7f4503480e7` | 35090 | `manuscript/figures/figure6_realdata_summary.png` |
| `d9a7ebf763bbd0fd920a3a9cacfc734d4c0eec9d7fd8bfeec592c1d12262762d` | 179191 | `manuscript/figures/figure7_group_block_average.png` |
| `776d5d7e75aa1bc27166045d751a76f7b4b15ab624268196b9590b59deee9a81` | 1978296 | `manuscript/main.pdf` |
| `8714e869f05e6ef2adf7baac55dcd2bf4a2bc8cc6f9e4bdce766606c16600537` | 174672 | `manuscript/main.tex` |
| `fd9a05f79869e6083d4279fdc222b89085d26aca07393fed95aa46fc5bca0b31` | 4120 | `reproduce_all.sh` |
| `01ce4001ad0fb9b37482eaadc6af393c6e0c5a5262b58533fe34c53802d98946` | 4112 | `requirements-lock.txt` |
| `90a1b0c629cceb68fb0900439a672b99a931e6f3ecba9ae2c8da9938a02eb9b0` | 1962 | `requirements.txt` |
| `2e06607221e8e81db393b5eb73f915f6c69e9ef95ba2844ab8f60441472fa0db` | 2701 | `results/fcortex_production.csv` |
| `d76657a7b24b5500796bc072ed2477610eb2dfb81f89135727ed00c2212ccb00` | 29980 | `results/fcortex_production.json` |
| `8f8b2e1e510fd893593186085908a957d8e6274a0e84a488776e30e7c9db4cd0` | 3685 | `results/homogeneous_fluence_validation.json` |
| `4b3424f183861894d3c966302c4baeb977e9e1097c462480c54d0cb7580e52ba` | 770 | `results/legacy/README.md` |
| `53ad26f40b8ae6bad8936b2cbb5eee008e28836d3835c6347878881bb9c2bd7f` | 771 | `results/legacy/mc_uncertainty_2layer.csv` |
| `1f574475f248a4361fb7b4e1f162c8cc8d7728f40934813fb0e4fe7a9b21aca8` | 6389 | `results/legacy/mc_uncertainty_2layer_and_csf.json` |
| `9f6092021e8515dc3527316eb626d5f6de852ccb3495d38edc2baee80f992019` | 2193 | `results/multiseed_operating_regime.json` |
| `ef259487e99aaae82229912bca5e03900f3945b900fd53a4a41baf3b757546f1` | 21509 | `results/realdata_v2_summary.json` |
| `6d6835b91063aee9f551ab796b40fa2c014906be0d85c00a21cc83b646fceb11` | 4985 | `results/robustness_secondary.json` |
| `eb8aab717fe206fc1bf96e987ed5981f3fdf4093cc6ca704740b714e442601a9` | 395807 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.pdf` |
| `2bbca4aa2b28f75a09ba812aa7183b533ccc494ba20f860ea9b85e488f5cf87a` | 108761 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.tex` |
| `928e138ed7d676fca0a6ecdcd61ec4a4c3556479b0b27865a617d350b2db00df` | 152680 | `supplementary/fnirs_kappa_beginner_notebook.ipynb` |

## Reproduce everything

Run `./reproduce_all.sh` (or `FAST=1 ./reproduce_all.sh` for a smoke run). Canonical commands:

- Forward model: `python code/mc_production.py -N 2000000 --batches 16 --thin_csf_N 2000000 --out fcortex_production --git_commit <sha> --analysis_round <label>` → copy JSON/CSV to `results/`.
- Robustness (Fig 3 source): `python code/mc_robustness_sweeps.py -N 500000 --seed 1 --git_commit <sha> --analysis_round <label> --out results/robustness_secondary.json`.
- Fluence validation: `python code/validate_homogeneous_fluence.py --out results/homogeneous_fluence_validation.json`.

`fcortex_source.py` requires schema 2.0, requires and re-verifies the payload SHA-256, and resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`).
