# Appendix C — Jang exercise/panic validation

## Question

Does accumulated action history or learned action-specific safety predict exercise, and does it improve next-day panic forecasting?

## Design and result

The public dataset contains daily smartphone and wearable measures from 43 outpatients with mood/anxiety-related diagnoses; 41 contribute complete consecutive-day rows. The fixed primary mapping treats exercise as the action and absence of panic on the following consecutive day as the favourable outcome. Prediction holds out each participant in turn.

Under the weak prior, the mechanism-informed exercise model improves on previous-action persistence (participant-level mean log-loss difference −0.0328, 95% bootstrap interval −0.0616 to −0.0075; better for 28/41 participants). A simpler action-history model is nevertheless best overall (log loss 0.243; Brier score 0.072); its advantage over the mechanism model is 0.020 log-loss units, although the participant-bootstrap interval includes zero (−0.010 to 0.054). The exact no-intercept specification is second worst (log loss 0.343), only slightly better than the intercept (0.354). As in the Fisher et al. dataset and analysis, overall performance is driven by persistence: the action-history model's change-only log loss is 1.418, versus 1.181 for the intercept.

For next-day panic, all 11 prespecified models are reported. Calibrated personal panic history performs best (log loss 0.107; Brier score 0.025), and adding the action-specific mechanism slightly worsens prediction. The intercept's below-chance leave-one-participant-out ranking statistics arise because each participant is predicted from the event rate among all other participants; these cells are not interpreted as ordinary in-sample intercept discrimination.

## Canonical files

- Time series: `../../data/raw/jang_openesm/0017_jang_ts.tsv`
- Static data: `../../data/raw/jang_openesm/0017_jang_static.tsv`
- Codebook: `../../data/raw/jang_openesm/0017_jang_codebook.pdf`
- Provenance: `../../provenance/zenodo_17347470_jang.json`
- Analysis: `../../run_validation.R`
- Action results: `../../results/jang_action_models.csv`
- Next-day panic results: `../../results/jang_event_models.csv`

## Interpretation limit

Exercise is uncommon and panic events are sparse in this dataset. The empirical action-conditional safety frequencies are only proxies for, not direct measurements of, the theoretical likelihood discriminabilities. The action model supports path-dependent behaviour, but personal event history is the stronger clinical forecasting signal. (Needless to say, observational prediction does not establish that exercise causally updates beliefs or reduces panic.)

## Source

Jang, S., et al. (2024). *A digital phenotyping dataset for impending panic symptoms: A prospective longitudinal study*. Scientific Data, **11**: 1264. https://doi.org/10.1038/s41597-024-04147-6
