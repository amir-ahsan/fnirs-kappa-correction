# Legacy results — NOT USED FOR THE MANUSCRIPT

These files are outputs of the superseded `mc_uncertainty.py` driver (an earlier,
lower-photon-count uncertainty workflow). They are retained only for historical
reference and are **not** the source of any manuscript number.

The authoritative forward-model artifacts are:

- `../fcortex_production.json` / `.csv` — produced by `code/mc_production.py`
  (schema 2.0, launch-defined paired batches, payload hash). Single source of
  truth for every f_cortex and CSF ratio in the manuscript.
- `../robustness_secondary.json` — produced by `code/mc_robustness_sweeps.py`
  (the secondary thickness / optical-property sweeps and Figure 3 source).

Do not use the files in this folder to reproduce the manuscript tables.
