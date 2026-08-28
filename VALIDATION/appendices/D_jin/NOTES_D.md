# Appendix D -- Jin sparse symptom trajectories

## Question

Can a reduced nonlinear transition map improve held-out prediction in the nine released weekly psychosis, depression, and mania histories, and (separately) what dynamical regimes are implied by the supplied Lorenz-model fits?

In other words, we are looking to distinguish forecasting evidence from dynamical interpretation.

## Design and result

The released subset contains nine selected patients, each represented by 209 reconstructed weeks across three binary symptom dimensions (psychosis, depression, and mania). Leave-one-patient-out comparisons include an intercept, first-order Markov persistence, linear symptom burden, a two-state hidden Markov model, and our active-inference-derived nonlinear map. The supplied Lorenz fits (from the original paper) are audited separately because they use a different in-sample inversion and observation model.

Only 34 of 5,401 finite adjacent same-symptom pairs change. The forecasting table contains 5,414 finite targets because models can issue an intercept or fallback forecast when the immediately preceding value of that symptom is missing. First-order Markov persistence dominates overall prediction (log loss 0.040; Brier score 0.007), whereas the nonlinear map scores 0.390 and 0.123. Both models are highly surprised by the rare changes (transition-only log loss 4.607 and 4.630, respectively, versus 0.756 for the intercept). Nonlinear coupling exceeds one in every fold, but only two of nine fitted fields lie within the bistable wedge. The supplied baseline Lorenz Rayleigh parameters remain below the classical stability threshold, with time-varying exceedance only for two patients.

## Canonical files

- Released data: `../../data/jin_2023/Data_for_Model.mat`
- Supplied Lorenz fits: `../../data/jin_2023/GCM_nosology.mat`
- Analysis: `../../../Python_code/clinical_trajectory_comparison.py`
- Explanatory notebook: `../../../Notebooks/Notebook_04_ClinicalTrajectories.md`
- Model comparison: `../../../outputs/clinical_trajectories/model_comparison.csv`
- Subject metrics: `../../../outputs/clinical_trajectories/subject_model_metrics.csv`
- Lorenz audit: `../../../outputs/clinical_trajectories/lorenz_fit_audit.csv`
- Figure: `../../../outputs/clinical_trajectories/clinical_trajectory_stress_test.png`

## Interpretation limit

The subset is selected, small, and dominated by unchanged weeks. Symptom burden is not engagement or withdrawal; fitted nonlinear parameters are not estimates of the manuscript's psychological controls. These data supply a useful negative boundary test, not validation of bistability.

## Source

Jin, J., Zeidman, P., Friston, K. J., & Kotov, R. (2023). *Inferring trajectories of psychotic disorders using dynamic causal modeling*. Computational Psychiatry, **7**(1), 60–75. https://doi.org/10.5334/cpsy.94
