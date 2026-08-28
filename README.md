# Action-Dependent Learning and Bistable Psychopathology Dynamics

This repository is the Supplementary Online Material (SOM) for a theory-led study of abrupt and path-dependent psychological change. The model closes a learning–action loop: behaviour determines which observations are sampled, sampled observations update policy-contingent beliefs, and those beliefs change subsequent behaviour.

The central formal quantity is the experience-dependent coupling function $\mathcal G(\tau;p)$. It rises and then falls as experience accumulates, so multistability is possible **only over an intermediate learning window**. The finite-agent analyses distinguish this quasistatic fixed-point structure from the slower stochastic trajectory that moves through it.

The evidence is ordered as follows:

1. **Formal and simulation results.** First, we show that the mean-field reduction permits folds and a local cusp normal form. The recovery protocol is evaluated at an intervention state with two stable branches and one unstable branch; all conditions are audited against the local critical boundary.
2. **Primary empirical trajectory.** In the Wichers dataset (N-of-1 observation series), a blind single-boundary search selects the same weekly occasion as the published analysis (day 127). The step is strongly supported against a constant mean but only modestly favoured over a delayed-discontinuation trend after penalising the breakpoint search.
3. **Secondary intensive analyses.** The Fisher dataset and the Jang dataset support behavioural-history prediction, but a simple action-history model generally outperforms the mechanism-informed model and prediction is weakest at action changes.
4. **Secondary sparse-trajectory analysis.** The released Jin-dataset histories strongly favour persistence overall, while all persistence-based models perform poorly at rare changes.

These results establish formal _plausibility_, empirical _compatibility_, and _partial support_ for path dependence. They do not uniquely identify the proposed action-contingent learning mechanism in clinical data.

## Start here

- [`SOM.md`](SOM.md): evidence map, appendix crosswalk, outputs, and claim limits.
- [`Notebooks/Notebook_01_ActiveInference.md`](Notebooks/Notebook_01_ActiveInference.md): formal derivation and the role of $\mathcal G(\tau;p)$.
- [`Notebooks/Notebook_02_MentalHealth.md`](Notebooks/Notebook_02_MentalHealth.md): clinical translation and the structurally valid recovery protocol.
- [`VALIDATION/README.md`](VALIDATION/README.md): public datasets, predictive design, and empirical limitations.
- [`VALIDATION/appendices/`](VALIDATION/appendices/): Appendix A–D routes matching the standalone appendices document.

## Repository structure

- `Python_code/`: formal utilities, onset simulation, bistable recovery analysis, conceptual figure, and Jin trajectory comparison.
- `Notebooks/`: derivation, clinical interpretation, and sparse-trajectory documentation.
- `VALIDATION/`: public intensive datasets, provenance, R analyses, results, and Appendix A–D notes.
- `outputs/manuscript_figures/`: canonical Figures 1–6 in PNG and PDF formats.
- `outputs/tables/`: onset and recovery audit tables.
- `outputs/clinical_trajectories/`: Jin model comparison, parameter audits, and figure.

## Complete reproduction

Python 3.10 or later, R, and the Python packages in `requirements.txt` are required. The distributed outputs were verified with Python 3.13.9 and R 4.4.1. The R validation uses base and recommended R packages only.

```bash
python3 -m pip install -r requirements.txt
python3 main.py
```

`main.py` regenerates the conceptual figure, onset audit, bistable recovery analyses, Jin comparison, Fisher/Jang/Wichers validation tables, and all six manuscript figures.

Individual routes are also available:

```bash
python3 Python_code/sim_psychopathology.py
python3 Python_code/recovery_bistable.py
python3 Python_code/clinical_trajectory_comparison.py
Rscript VALIDATION/run_validation.R
Rscript VALIDATION/generate_wichers_figure.R
```

## Important to bear in mind: Interpretation boundaries on empirical results

- The Wichers change-point endpoint contains 28 weekly depression scores. The surrounding approximately 1,476 momentary records establish study intensity, not 1,476 independent depression outcomes.
- Day 127 is selected in 86.9% of a 10,000-replicate residual bootstrap; the 95% highest-frequency set is {120, 127}.
- A maximum-F test rejects a constant mean, but the residual test after linear detrending does not distinguish a step from a monotone rise.
- The Fisher and Jang proxies are action-conditional outcome frequencies, not direct measurements of the theoretical likelihood discriminabilities.
- In the secondary datasets, strong aggregate prediction mainly reflects persistence and degrades substantially at changes.
- No included dataset randomised actions or directly manipulated the proposed learning loop.

## Version 2

- This is v2 of the analysis and codebase; the original work used the Millennium Cohort Study (MCS) and the Adolescent Brain Cognitive Development (ABCD) study as empirical datasets. The focus of those analyses was on the initially derived _balance index_. Those results were consistent with the corresponding corollary of the theory However, the empirical operationalization relied on proxy variables and a log-ratio quantity rather than directly examining the learning–action mechanism proposed in the main manuscript.
- Therefore, in this version we use a series of other, more suitable, empirical observations that are more directly aligned with quantities and predictions arising from the proposed mechanism. The analyses are thus refocused on empirical tests that remain closer to the mechanism described in the main manuscript and accompanying Notebooks in this SOM.
- Finally, the codebase/outputs and accompanying summaries in this version were additionally re-checked and stress-tested with assistance from OpenAI Codex using GPT-5.6 Sol (July 2026 Codex deployment; accessed August 2026). This review included checks for implementation errors, internal inconsistencies, and sensitivity of the reported conclusions to analysis choices. All analyses and resulting claims were subsequently inspected and verified.