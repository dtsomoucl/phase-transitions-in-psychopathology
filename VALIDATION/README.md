# Public trajectory validation

This folder contains the empirical programme re-analysing data from public sources. The analyses are ordered by what each dataset can identify; they are not treated as interchangeable replications.

## Primary result: Wichers N-of-1 trajectory

This public study followed one person with recurrent major depression for 239 days during venlafaxine discontinuation, with up to ten momentary assessments per day. The change-point endpoint is the 28-score weekly SCL-90-R depression series.

- Without using the originally reported boundary as an input, the single-boundary search selects experiment day 127, the same weekly occasion reported in the original analysis.
- A 10,000-replicate residual bootstrap selects day 127 in 86.9% of samples; the 95% highest-frequency set is {120, 127}.
- A maximum-F permutation test rejects a constant-mean null (*p* = .00010).
- After linear detrending, the maximum-F permutation result is not significant (*p* = .704).
- All model comparisons use the same 27 observations. Charging the step one additional breakpoint parameter gives adjusted BIC 11.39, compared with 13.23 for delayed discontinuation and 14.84 for linear time.

The key result is that the weekly series supports a **discontinuity relative to a constant mean** and the search recovers **the same weekly occasion**. However, clearly, the data do not cleanly distinguish a step from a monotone rise and do not identify the proposed learning mechanism.

## Secondary result: Fisher next-action prediction

Forty participants completed approximately four assessments daily for about 30 days. Engagement is defined as the mean of the two avoidance items (activity and people) below 50; favourable outcome is contentment at or above 50.

The mechanism-informed model is slightly better than previous-action persistence under the weak prior, but its bootstrap interval includes zero. The simple action-history model is best and outperforms the mechanism model (mean log-loss difference 0.058, 95% participant-bootstrap interval 0.035 to 0.086). Performance is driven by persistence and deteriorates at changes.

## Secondary result: Jang exercise and panic

Forty-one participants contribute complete consecutive-day rows. Exercise is the action and no panic on the following day is the favourable outcome.

The mechanism-informed action model improves on previous-action persistence under the weak prior, but the action-history model remains best. Its advantage over the mechanism model is small and uncertain (0.020, 95% interval −0.010 to 0.054). For next-day panic, calibrated personal panic history performs best; adding the action-specific mechanism slightly worsens prediction. All 11 prespecified panic models are reported.

## Construct mapping limit

For Fisher and Jang, the quantities used as learned discriminabilities are action-conditional outcome frequencies. They are not direct measurements of the likelihood discriminabilities whose entropies generate the ambiguity term in the theoretical model. These analyses test behavioural history dependence and a reduced predictive mapping, not full parameter recovery.

## Reproduction

From the repository root:

```bash
Rscript VALIDATION/run_validation.R
Rscript VALIDATION/generate_wichers_figure.R
```

The scripts use base and recommended R packages. They write to `VALIDATION/results/` and `outputs/manuscript_figures/`.

To print results without replacing the distributed tables:

```bash
VALIDATION_WRITE=0 Rscript VALIDATION/run_validation.R
```

## Generated outputs

- `results/fisher_action_models.csv`
- `results/jang_action_models.csv`
- `results/jang_event_models.csv`
- `results/wichers_transition_models.csv`
- `results/wichers_transition_summary.csv`
- `results/wichers_change_point_inference.csv`
- `results/wichers_boundary_bootstrap.csv`
- `results/wichers_daily_density.csv`
- `results/wichers_weekly_series.csv`
- `results/wichers_change_point_profile.csv`
- `../outputs/manuscript_figures/Fig_6_Wichers.png`
- `../outputs/manuscript_figures/Fig_6_Wichers.pdf`

## Integrity records

- Fisher time-series MD5: `04dbe56452b133242626a2ec826d0ad9`
- Fisher codebook MD5: `ba0d5f406435409da711d68405260a39`
- Jang time-series MD5: `6050a59989968cb4606480f6b8ee02c9`
- Jang static-data MD5: `e7036c95cd8db7fae3fd50751f2acf7b`
- Jang codebook MD5: `4a0a4fe784b2d1877353c74e84679131`
- Wichers archive MD5: `d17d37f0df5d715457f12b38b6e8b186`
- Wichers archive SHA-256: `465ae9862f6d8a1d10ff322dfc1a88a87b422988501555d9cbedd85fdb535c5f`

Machine-readable source metadata are retained under `provenance/`.
