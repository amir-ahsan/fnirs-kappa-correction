# Release Manifest

Production data provenance (from `results/fcortex_production.json`):

- `schema_version`: 2.0
- `git_commit` (data, real SHA): 32687ed  |  `git_commit_full`: 32687ede2ff6786b0be069a84e80266e8ee5f321  |  `git_dirty`: False
- `analysis_round`: round9
- `generated_utc`: 2026-08-07T10:01:24Z
- `data_sha256` (numeric payload): `91cfa828227df0652fb023a9ec5a005424103bcc5841c162fe9cfc58d2ac0331`
- SDS grid (mm): [25.0, 30.0, 33.0, 35.0, 38.0, 40.0, 42.0]  (spans the in-vivo channel range 33.4-40.9 mm)
- `N_per_config`: 2000000, `n_batches`: 16, `seed`: 1, `N_thin_csf`: 2000000, `g`: 0.9, `L_max`: 1200.0, `z_max`: 150.0
- `python_version`: 3.11.15  |  platform: Linux-6.18.5-fc-v18-x86_64-with-glibc2.39
- packages: numpy 2.4.4, scipy 1.17.1, pandas 3.0.2, matplotlib 3.10.9, mne 1.12.1, mne-bids 0.19.0, h5py 3.16.0

All artifacts share the uniform provenance schema written by `code/provenance.py` (git_commit / git_commit_full / git_dirty are auto-detected from HEAD; a human release name is stored separately as `release_label`). Each artifact records the commit it was generated from in its own `_meta`. The production forward model and the secondary MC artifacts are frozen at `round9`/`32687ed`; `results/realdata_v2_summary.json` was regenerated later (`round13`) to record the truthful per-run dataset acquisition provenance and content (tree) hash, using the unchanged `round9` production (its `input_hashes.fcortex_production_sha256` still points to the same production payload). The numerical results are unchanged.

## SHA-256 of tracked files

| SHA-256 | Size (bytes) | File |
|---------|-------------:|------|
| `0bcc9c77d89fcdbf3b849ad8832b283a4519c4c5c6bfaf08c48c97723d3abc5f` | 484 | `.gitignore` |
| `1c941fc2d0be0abea275e87cb5a86649e7a5f227f1cf095da98812bcc504d4c0` | 1045 | `CITATION.cff` |
| `9d020d65a2756c2fc42cbbfdd7d88151cefcc1407cd1dad6814e73f5807f7d05` | 1083 | `LICENSE` |
| `0050fc32fc270930da596b499b58e5bda4da0f3a2ce2827e95938a7e52d99672` | 953 | `LICENSE-manuscript.md` |
| `c1c6a2496cb08fc785c73c0ebe39a64268971dafd832a918c1fb16fd2288ad4c` | 17204 | `README.md` |
| `4724aedbf26e88b84b434eed8b67a2535920df3b960ff69cc867cc68fd080a1e` | 106430 | `REVISION_RESPONSE.md` |
| `1f9705d3eba85093b106f62b46e08bc7e10c9bf44c95bbc396bf85a1f39fa75b` | 6985 | `code/fcortex_source.py` |
| `7584b084574a784641b943ae7acd5b17a66a1248d376c45b6d718115c461d9a5` | 20815 | `code/fnirs_invivo_demo.py` |
| `546eb0b81f3dc9a652063a31da30f06450e0c9f58c6cddaf3b6255cb39c8b4a6` | 14625 | `code/fnirs_kappa_group_analysis.py` |
| `76c45fe8c45a1f30a9a8bcba687b9500fc872054733b5bbf50d72177b832e92b` | 46518 | `code/fnirs_kappa_realdata_analysis.py` |
| `b9fd8616c397a79e73fab2f2dba19e96a46e1adc23565b43d1520e9a8c88c435` | 29575 | `code/fnirs_kappa_realdata_v2.py` |
| `4c3a17184641ba788a2d0fac324aae10f90d05e0ddd32049bc27ba0f3b2ddc07` | 108029 | `code/fnirs_kappa_synthetic_validation.py` |
| `42eae86326cf757b6f1a9ac956db772671de75efa1b2f53870621e083bf20a2c` | 5013 | `code/make_release_manifest.py` |
| `35e1fe777e0dc3657b5ba313579674cd6eb99d418b00614b5f1f583f16c919be` | 9364 | `code/mc_2layer.py` |
| `4544749aabea9c2834e72610cdeaba2195b5c11594a76ba959407b203b484dfe` | 9408 | `code/mc_csf.py` |
| `30b0f092a53c4b3b5e0768c4f2dd56cd2353a150822d91aefcb8ab13223aefec` | 23234 | `code/mc_production.py` |
| `e2607c7c4efb882c968843c57373d5b32b1f8c3c5179ff3c5056b7a2a9568868` | 9855 | `code/mc_robustness_sweeps.py` |
| `4be23d4f324bab85fbe7ec831930bfd8bfdc3926cb683d2143531eff8bd90d04` | 7603 | `code/mc_uncertainty.py` |
| `8b18f5868f82318fc12e5be67ceb470b6b55710d213dffcf5aac3d3a45b67f9b` | 5363 | `code/provenance.py` |
| `f5c75ecf36755d119c93ddf523555c372c4cce46548013c219935a3f6d98ae7c` | 8290 | `code/validate_homogeneous_fluence.py` |
| `38f2ef56e46caeb5726141d7714bbc49ab2add6f744fc87be8cb888f06377990` | 568604 | `manuscript/figures/figure1_timeseries.png` |
| `20d522737804c3f4cd02741090503e53bdc9a5947f6fb9cf8a61c4a1f28a2aea` | 226267 | `manuscript/figures/figure2_summary.png` |
| `95b6b39f7350a1bd8eedcb235ddc80d51f703116f391d75614bec13a57689ad6` | 533977 | `manuscript/figures/figure3_robustness.png` |
| `a56ecf9cdabee52dbcfb1da8ec70c08dde1b02fd8085bae6489bc5fbced1c26a` | 64405 | `manuscript/figures/figure4_realdata_timeseries.png` |
| `1830309e44230290c6c25ebacb6ee6ee36f3b1447c3600317b4686337b9852a5` | 96158 | `manuscript/figures/figure5_realdata_hrf.png` |
| `8eee074082111c72b5f38b02b3668209c37a488afbdebe51276387c7c6fcfc22` | 35028 | `manuscript/figures/figure6_realdata_summary.png` |
| `e86640dc8f9999189a8c9570c21f7d1acd8b94c4c6ff6bdce3767351e0f3f73e` | 178654 | `manuscript/figures/figure7_group_block_average.png` |
| `2602e27d603d95f77729174471362a691c7592d45ada2bc099b35eeeabde0c5d` | 1997120 | `manuscript/main.pdf` |
| `d653f97adb08d81bd0f9cd6c69691b4203eeafe743fa471094ddc35f37842e58` | 189811 | `manuscript/main.tex` |
| `6241a863070d117b80ac2d28df38021fd1ebd2127d21f5e41b973667b8d79c68` | 6364 | `reproduce_all.sh` |
| `f50e014ad9afaa7cd8a0d1c79b91d516570ac87e72ce4c8d052f862ba7d793a1` | 2687 | `requirements-lock.txt` |
| `90a1b0c629cceb68fb0900439a672b99a931e6f3ecba9ae2c8da9938a02eb9b0` | 1962 | `requirements.txt` |
| `2e06607221e8e81db393b5eb73f915f6c69e9ef95ba2844ab8f60441472fa0db` | 2701 | `results/fcortex_production.csv` |
| `ec21bb6e18b0355b5df86f018e189e8f3d0cc7d97efe02af1fea03f2a0a2f1ad` | 30610 | `results/fcortex_production.json` |
| `9b59a3e4e569e83fc6af26912d8aeb6a99c9963565f7042b01b4289fc018d09e` | 3664 | `results/homogeneous_fluence_validation.json` |
| `4b3424f183861894d3c966302c4baeb977e9e1097c462480c54d0cb7580e52ba` | 770 | `results/legacy/README.md` |
| `53ad26f40b8ae6bad8936b2cbb5eee008e28836d3835c6347878881bb9c2bd7f` | 771 | `results/legacy/mc_uncertainty_2layer.csv` |
| `1f574475f248a4361fb7b4e1f162c8cc8d7728f40934813fb0e4fe7a9b21aca8` | 6389 | `results/legacy/mc_uncertainty_2layer_and_csf.json` |
| `a2a44c3fcb81d2dfeab110d62196502bd4b8c02907feafdfd3716d5c84f1f3e4` | 2149 | `results/multiseed_operating_regime.json` |
| `0bf6e42d58260bde0e0e36d17049706528aa9c55f9e0f6fd035e53268798b523` | 33348 | `results/realdata_v2_summary.json` |
| `372b515b10c2ee116f713434d0553196d32d4e9f377b1b34f36def874fee018b` | 4964 | `results/robustness_secondary.json` |
| `86f8199e1f7bdf8df289eabefca96eff14301d0551e45355806b62eb826487d5` | 440288 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.pdf` |
| `9ea7431ca3a7475566b172b4f50d8ba7edaf470b73e9cb4d8d13b5d6ff515e80` | 117891 | `supplementary/fNIRS_Kappa_Pedagogical_Guide.tex` |
| `133034436a00abd23002f0fad223eb96966c9705664cd85b0fc0915efffcd7e2` | 150124 | `supplementary/fnirs_kappa_beginner_notebook.ipynb` |

## Reproduce everything

Run `./reproduce_all.sh` (or `FAST=1 ./reproduce_all.sh` for a smoke run; `RELEASE=1 ./reproduce_all.sh` makes the in-vivo step fatal and regenerates this manifest automatically). Canonical commands:

- Forward model: `python code/mc_production.py -N 2000000 --batches 16 --thin_csf_N 2000000 --out fcortex_production --analysis_round <label>` → copy JSON/CSV to `results/`.
- Robustness (Fig 3 source): `python code/mc_robustness_sweeps.py -N 500000 --seed 1 --analysis_round <label> --out results/robustness_secondary.json`.
- Fluence validation: `python code/validate_homogeneous_fluence.py --out results/homogeneous_fluence_validation.json`.

`fcortex_source.py` requires schema 2.0, requires and re-verifies the payload SHA-256, and resolves the tracked `results/` copy first (override with `FCORTEX_PRODUCTION_JSON`).
