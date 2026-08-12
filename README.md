# Active-Inference Phase Transitions in Psychopathology

Code and analysis materials for a project that combines:

- theoretical and simulation work on phase transitions in engagement versus withdrawal (anchored in adolescent mental health)
- empirical longitudinal analyses in the Millennium Cohort Study (MCS) and Adolescent Brain Cognitive Development Study (ABCD)
- all code needed for figures and outputs used in the study

## Repository structure

- `Notebooks/`
  - Conceptual notebooks that introduce the core derivations and model logic.
  - `Notebook_01_ActiveInference.md`: main active-inference and phase-transition derivations
  - `Notebook_02_MentalHealth.md`: mental-health interpretation and clinical framing
  - `Notebook_03_BalanceIndex.md`: exploratory ratio-composite derivation and scale assumptions

- `Python_code/`
  - Main simulation code for the theoretical and computational results.
  - `sim_psychopathology.py`: onset simulations and core transition figures
  - `step2_hysteresis.py`: recovery-boundary and field-dominance analyses
  - `step3_orthogonal.py`: orthogonal intervention analyses
  - `step4_decay.py`: symmetric-decay analyses
  - `step5_asymmetric_memory.py`: asymmetric-retention extension
  - `step6_verified_withdrawal.py`: fixed-schedule state diagnostic and paired state-anchored recovery sensitivity
  - `export_manuscript_figures.py`: exports manuscript versions of simulation figures
  - `core_functions.py`, `psychopathology_regimes.py`: shared model functions and parameter regimes
  - `empirical_bridge.py`: calibrated synthetic demonstration linking simulation concepts to empirical constructs; not independent evidence (i.e. computational phenotyping remains to be done in future work)

- `R/`
  - Main empirical analysis code for the cohort-study results.
  - `abcd_pipeline.R`: ABCD data-processing and model pipeline
  - `mcs_pipeline.R`: MCS data-processing and model pipeline
  - `models.R`: shared model fitting, summary, and plotting functions
  - `balance_index.R`: exploratory ratio-composite construction and evaluation
  - `reliability.R`: scale reliability checks
  - `inventory.R`, `utils.R`: inventorying and shared utilities

- `scripts/`
  - `run_pipeline.R`: top-level R entry point for the empirical workflow
  - `regenerate_figure6.R`: regenerates the cross-cohort coefficient figure from aggregate tables
  - `rerun_abcd_fiml.R`: optional targeted rerun of the ABCD FIML sensitivity from an authorised analysis cache

- Repository root
  - `main.py`: top-level Python entry point for the full simulation workflow
  - `SOM_Empirical_Results.Rmd`: main supplemental empirical-results document

- `config/`
  - Project configuration, including data and output locations

---

## Software requirements

The simulation workflow uses Python 3.10 or later. Install its dependencies with:

```bash
python3 -m pip install -r requirements.txt
```

The empirical workflow and SOM use R 4.3 or later with the following packages: `dplyr`, `fs`, `ggplot2`, `glue`, `gt`, `haven`, `jsonlite`, `kableExtra`, `knitr`, `lavaan`, `lme4`, `lmerTest`, `mice`, `mitools`, `psych`, `purrr`, `readr`, `readxl`, `rmarkdown`, `sandwich`, `scales`, `stringr`, `survey`, `tibble`, `tidyr`, and `yaml`. Rendering the SOM also requires Pandoc.

## Suggested order of work (for reproducibility)

### 1. Read the notebooks first

For orientation, start with the notebooks in `Notebooks/`, especially:

1. `Notebook_01_ActiveInference.md`
2. `Notebook_02_MentalHealth.md`
3. `Notebook_03_BalanceIndex.md`

These notebooks explain the derivations, the mapping to phase-transition mathematics, and how the empirical proxies were chosen.

(Related/background work can also be accessed in a different repo, which focuses on developmental psyschology context: https://github.com/dtsomoucl/phase-transitions-in-active-inference)

### 2. Run the Python simulations

From the repository root:

```bash
python3 main.py
```

This generates the main simulation outputs and figures under `Python_code/Figs_psychopathology/`.

To export manuscript versions of the simulation figures:

```bash
python3 Python_code/export_manuscript_figures.py
```

The exporter preserves the manuscript numbering: conceptual overview (Fig. 1), onset (Fig. 2), intervention comparison (Fig. 3), recovery boundary (Fig. 4), and retention analysis (Fig. 5). Figure 6 is regenerated separately from aggregate empirical coefficient tables with `Rscript scripts/regenerate_figure6.R`.

### 3. Run the empirical pipeline

From the repository root:

```bash
Rscript scripts/run_pipeline.R
```

This runs the MCS and ABCD empirical workflows, assuming authorised users have separately obtained the licensed data from the UK Data Service and NIH and configured the protected input directories. Licensed data are not included in this repository. The pipeline writes outputs to:

- `outputs/derived/`
- `outputs/tables/`
- `outputs/figures/`

### 4. Render the supplemental results

The repository includes the non-disclosive aggregate tables and figures required to render the SOM without licensed data or participant-level analysis caches:

```bash
Rscript -e "rmarkdown::render('SOM_Empirical_Results.Rmd')"
```

For full reproducibility, one would need to get the correct license agreement and access/use the MCS and ABCD datasets (the releases used here, for this research project, or the latest ones but this would be subject to minor differences potentially).

## Outputs

Main generated outputs include:

- simulation figures in `Python_code/Figs_psychopathology/`
- empirical tables and figures in `outputs/tables/` and `outputs/figures/`
- manuscript-specific figures in `Manuscript/`
- state-anchored sensitivity tables in `outputs/tables/fixed_schedule_intervention_state.csv` and `outputs/tables/verified_withdrawal_*.csv`
- the state-anchored sensitivity figure `Python_code/Figs_psychopathology/fig_S6_verified_withdrawal.*`

## Notes

- The Python and R codebases are complementary: the Python code establishes the theoretical and simulation results, while the R code evaluates a narrower observational ordering of predictor families in longitudinal cohort data. The cohort regressions do not estimate the model's latent parameters or identify its transition mechanism.
- The notebooks are the best place to understand the model assumptions and the logic linking the formal theory to the empirical analyses (and they also contain some exploratory material).
- In legacy function arguments (which were used in the initial version of the project), `T_ill` denotes elapsed trials on the fixed adverse schedule. It is not verified illness duration or withdrawal-residence time. The main field-versus-precision ordering survives the state-anchored analysis, whereas the asymmetric-memory duration pattern does not generalise monotonically to additional trials after verified withdrawal.
- Rendering `SOM_Empirical_Results.Rmd` uses only the distributed aggregate tables and figures. It does not load licensed cohort files or participant-level analysis caches.


D. I. Tsomokos, 12/Aug/2026