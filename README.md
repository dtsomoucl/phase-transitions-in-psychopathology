# Active-Inference Phase Transitions in Psychopathology

Code and analysis materials for a project that combines:

- theoretical and simulation work on phase transitions in engagement versus withdrawal (anchored in adolescent mental health)
- empirical longitudinal analyses in the Millennium Cohort Study (MCS) and Adolescent Brain Cognitive Development Study (ABCD)
- all figures and outputs used in the study

## Repository structure

- `Notebooks/`
  - Conceptual notebooks that introduce the core derivations and model logic.
  - `Notebook_01_ActiveInference.md`: main active-inference and phase-transition derivations
  - `Notebook_02_MentalHealth.md`: mental-health interpretation and clinical framing
  - `Notebook_03_BalanceIndex.md`: balance-index derivation and interpretation

- `Python_code/`
  - Main simulation code for the theoretical and computational results.
  - `main.py`: top-level Python entry point for the simulation workflow
  - `sim_psychopathology.py`: onset simulations and core transition figures
  - `step2_hysteresis.py`: recovery-boundary and field-dominance analyses
  - `step3_orthogonal.py`: orthogonal intervention analyses
  - `step4_decay.py`: symmetric-decay analyses
  - `step5_asymmetric_memory.py`: asymmetric-retention extension
  - `export_manuscript_figures.py`: exports manuscript versions of simulation figures
  - `core_functions.py`, `psychopathology_regimes.py`: shared model functions and parameter regimes
  - `empirical_bridge.py`: bridge code linking simulation concepts to empirical constructs

- `R/`
  - Main empirical analysis code for the cohort-study results.
  - `abcd_pipeline.R`: ABCD data-processing and model pipeline
  - `mcs_pipeline.R`: MCS data-processing and model pipeline
  - `models.R`: shared model fitting, summary, and plotting functions
  - `balance_index.R`: balance-index construction and evaluation
  - `reliability.R`: scale reliability checks
  - `inventory.R`, `utils.R`: inventorying and shared utilities

- `scripts/`
  - `run_pipeline.R`: top-level R entry point for the empirical workflow

- `config/`
  - Project configuration, including data and output locations

---

## Typical workflow

### 1. Read the notebooks first

For orientation, start with the notebooks in `Notebooks/`, especially:

1. `Notebook_01_ActiveInference.md`
2. `Notebook_02_MentalHealth.md`
3. `Notebook_03_BalanceIndex.md`

These notebooks explain the derivations, the mapping to phase-transition mathematics, and how the empirical proxies were chosen.

### 2. Run the Python simulations

From the repository root:

```bash
python3 Python_code/main.py
```

This generates the main simulation outputs and figures under `Python_code/Figs_psychopathology/`.

To export manuscript versions of the simulation figures:

```bash
python3 Python_code/export_manuscript_figures.py
```

### 3. Run the empirical pipeline

From the repository root:

```bash
Rscript scripts/run_pipeline.R
```

This runs the MCS and ABCD empirical workflows (assuming you have obtained access to the data under license from the UK Data Service and the NIH, and have placed the data under folders data_MCS and data_ABCD), and writes outputs to:

- `outputs/derived/`
- `outputs/tables/`
- `outputs/figures/`

## Outputs

Main generated outputs include:

- simulation figures in `Python_code/Figs_psychopathology/`
- empirical tables and figures in `outputs/tables/` and `outputs/figures/`
- manuscript-specific figures in `Manuscript/`

## Notes

- The Python and R codebases are complementary: the Python code establishes the theoretical and simulation results, and the R code tests the corresponding empirical ordering predictions in longitudinal cohort data.
- The notebooks are the best place to understand the model assumptions and the logic linking the formal theory to the empirical analyses (and they also contain some exploratory material).
