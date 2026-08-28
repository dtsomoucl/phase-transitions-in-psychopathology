# Appendix A — Wichers N-of-1 abrupt transition

## Question

Can a transparent analysis of a well-characterised intensive clinical series recover the abrupt transition previously identified at day 127 (i.e. as published in the original article)?

## Design and result

The public study followed one 57-year-old man with recurrent major depression for 239 days while venlafaxine was tapered from 150 mg to zero, with up to ten momentary assessments per day. The present change-point endpoint is the 28-score weekly SCL-90-R depression series. A data-driven single-boundary search, constrained to leave at least four weekly observations on either side, selects the same weekly occasion as the original publication: experiment day 127.

All comparison models were fitted to the common 27-observation subset. Because the published boundary was derived from the same data, the step model was charged one additional breakpoint parameter. Its adjusted BIC was 11.39, compared with 13.23 for a delayed-discontinuation model, 14.84 for linear time, and 32.07 for contemporaneous linear medication dose. The step is therefore modestly (not decisively) favoured over the strongest gradual comparator. Mean weekly depression was 1.255 before and 1.981 after the boundary. Among observations at the same zero medication dose, the mean was 1.354 before versus 1.981 after the boundary (difference 0.627; 5 versus 12 weekly observations).

A 10,000-replicate residual bootstrap selected day 127 in 86.9% of samples; the 95% highest-frequency set was {120, 127}. A maximum-F permutation test rejected a constant-mean null (*p* = .00010), whereas the same test applied after linear detrending did not distinguish a step from residual fluctuation (*p* = .704).

## Canonical files

- Source archive: `../../data/raw/wichers_transition/ESMdata.zip`
- Provenance: `../../provenance/osf_j4fg8_files.json`
- Analysis: `../../run_validation.R`
- Figure generator: `../../generate_wichers_figure.R`
- Model comparison: `../../results/wichers_transition_models.csv`
- Summary: `../../results/wichers_transition_summary.csv`
- Daily assessment density: `../../results/wichers_daily_density.csv`
- Weekly endpoint series: `../../results/wichers_weekly_series.csv`
- Change-point profile: `../../results/wichers_change_point_profile.csv`
- Permutation tests: `../../results/wichers_change_point_inference.csv`
- Boundary bootstrap: `../../results/wichers_boundary_bootstrap.csv`
- Manuscript figure: `../../../outputs/manuscript_figures/Fig_6_Wichers.png`

## Interpretation limit

This is a reproducible reanalysis of the same public series, not an independent replication. The weekly endpoint is embedded in a dense experience-sampling design, but the change-point search does not use approximately 1,476 independent depression outcomes. The result supports a discontinuity relative to a constant mean and recovers the same weekly occasion, but it does not cleanly separate a step from a monotone rise and cannot identify the proposed action-contingent learning mechanism.

## Sources

- Kossakowski, J. J., Groot, P. C., Haslbeck, J. M. B., Borsboom, D., & Wichers, M. (2017). *Data from Critical Slowing Down as a Personalized Early Warning Signal for Depression*. Journal of Open Psychology Data, **5**(1). https://doi.org/10.5334/jopd.29
- Wichers, M., Groot, P. C., Psychosystems, ESM Group, & EWS Group. (2016). *Critical slowing down as a personalized early warning signal for depression*. Psychotherapy and Psychosomatics, **85**(2), 114–116. https://doi.org/10.1159/000441458
