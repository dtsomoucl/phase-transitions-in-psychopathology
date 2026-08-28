# Appendix B — Fisher avoidance/contentment validation

## Question

Does accumulated behavioural history, or outcome-contingent learning, predict subsequent engagement in an intensive depression/anxiety dataset?

## Design and result

Forty participants with generalised anxiety disorder (GAD), major depressive disorder (MDD), or both, completed ~4 assessments daily for roughly 30 days. The fixed primary mapping defines engagement as the mean of the two avoidance items (activity and people) below 50 and a favourable outcome as contentment at or above 50. Prediction holds out each participant in turn.

With weak Beta–Bernoulli priors, the mechanism-informed model is only slightly better than a previous-action model (participant-level mean log-loss difference −0.0056, 95% bootstrap interval −0.0185 to 0.0078; better for 22/40 participants). A simpler model using previous action plus cumulative action history predicts best (log loss 0.471; Brier score 0.152), and it outperforms the mechanism-informed model by 0.058 log-loss units (95% bootstrap interval 0.035 to 0.086). The mechanism model's advantage over previous action becomes favourable under stronger priors, but it remains worse than action history.

The exact expected-free-energy specification lacks an intercept and cannot reproduce the engagement base rate; it is therefore a specification diagnostic rather than a fair stand-alone mechanism test (log loss 0.696 versus 0.658 for the intercept). Across models, strong overall scores are driven by persistence: the best action-history model has change-only log loss 1.221, compared with 0.738 for the intercept.

## Canonical files

- Time series: `../../data/raw/fisher_openesm/0033_fisher_ts.tsv`
- Codebook: `../../data/raw/fisher_openesm/0033_fisher_codebook.docx`
- Provenance: `../../provenance/zenodo_17348039_fisher.json`
- Analysis: `../../run_validation.R`
- Results: `../../results/fisher_action_models.csv`

## Interpretation limit

Avoidance and contentment are measured in the same assessment window, so their temporal order is not established. Binary thresholds simplify continuous measures, and no action was randomised. The empirical quantities called learned discriminability are action-conditional outcome frequencies, not direct measurements of the likelihood discriminabilities whose entropies produce the theoretical ambiguity term. The result supports behavioural history dependence, not causal identification of outcome-contingent learning.

## Source

Fisher, A. J., Reeves, J. W., Lawyer, G., Medaglia, J. D., & Rubel, J. A. (2017). *Exploring the idiographic dynamics of mood and anxiety via network analysis*. Journal of Abnormal Psychology, **126**(8), 1044–1056. https://doi.org/10.1037/abn0000311
