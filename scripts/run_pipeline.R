### DT --> End-to-end runner for the local ABCD/MCS cohort pipeline.
source("R/utils.R")
source("R/inventory.R")
source("R/abcd_pipeline.R")
source("R/mcs_pipeline.R")
source("R/models.R")
source("R/reliability.R")

config <- load_project_config()

ensure_directories(unlist(config$outputs))

abcd_mapping <- abcd_variable_mapping(config)
mcs_mapping <- mcs_variable_mapping(config)

write_inventories(config, abcd_mapping, mcs_mapping)

abcd <- build_abcd_dataset(config)
mcs <- build_mcs_dataset(config)

safe_write_rds(abcd, fs::path(config$outputs$derived_dir, "abcd_analysis_objects.rds"))
safe_write_rds(mcs, fs::path(config$outputs$derived_dir, "mcs_analysis_objects.rds"))

abcd_desc <- describe_longitudinal_data(abcd$long, "participant_id", "session_id", "internalising_t", "ABCD")
mcs_desc <- describe_longitudinal_data(mcs$long, "MCSID", "wave", "distress_z", "MCS")

safe_write_csv(abcd_desc$wave_summary, fs::path(config$outputs$tables_dir, "abcd_wave_summary.csv"))
safe_write_csv(abcd_desc$missingness, fs::path(config$outputs$tables_dir, "abcd_missingness.csv"))
safe_write_csv(abcd_desc$attrition, fs::path(config$outputs$tables_dir, "abcd_attrition.csv"))

safe_write_csv(mcs_desc$wave_summary, fs::path(config$outputs$tables_dir, "mcs_wave_summary.csv"))
safe_write_csv(mcs_desc$missingness, fs::path(config$outputs$tables_dir, "mcs_missingness.csv"))
safe_write_csv(mcs_desc$attrition, fs::path(config$outputs$tables_dir, "mcs_attrition.csv"))

mcs_reliability <- build_mcs_composite_reliability(mcs)
abcd_reliability <- build_abcd_composite_reliability(abcd)

safe_write_csv(mcs_reliability, fs::path(config$outputs$tables_dir, "mcs_composite_reliability.csv"))
safe_write_csv(abcd_reliability, fs::path(config$outputs$tables_dir, "abcd_composite_reliability.csv"))

abcd_trajectory <- fit_trajectory_model(
  abcd$long %>% filter(!is.na(internalising_t), !is.na(wave_year), session_id <= config$analysis$abcd_followup_session),
  id_var = "participant_id",
  time_var = "wave_year",
  outcome_var = "internalising_t"
)

mcs_trajectory <- fit_trajectory_model(
  mcs$long %>% filter(!is.na(distress_z), !is.na(age_years)),
  id_var = "MCSID",
  time_var = "age_years",
  outcome_var = "distress_z"
)

safe_write_csv(tidy_lmer_summary(abcd_trajectory, "abcd_trajectory"), fs::path(config$outputs$tables_dir, "abcd_trajectory_model.csv"))
safe_write_csv(tidy_lmer_summary(mcs_trajectory, "mcs_trajectory"), fs::path(config$outputs$tables_dir, "mcs_trajectory_model.csv"))

abcd_results <- fit_abcd_models(abcd)
mcs_results <- fit_mcs_models(mcs)

safe_write_csv(abcd_results$primary_demographics, fs::path(config$outputs$tables_dir, "abcd_primary_demographics.csv"))
safe_write_csv(mcs_results$primary_demographics, fs::path(config$outputs$tables_dir, "mcs_primary_demographics.csv"))

safe_write_csv(abcd_results$metrics, fs::path(config$outputs$tables_dir, "abcd_model_metrics.csv"))
safe_write_csv(abcd_results$coefficients, fs::path(config$outputs$tables_dir, "abcd_model_coefficients.csv"))
safe_write_csv(abcd_results$moderation, fs::path(config$outputs$tables_dir, "abcd_context_moderation.csv"))
safe_write_csv(abcd_results$moderation_full, fs::path(config$outputs$tables_dir, "abcd_context_full.csv"))
safe_write_csv(abcd_results$coefficient_tests, fs::path(config$outputs$tables_dir, "abcd_coefficient_tests.csv"))
safe_write_csv(abcd_results$sensitivity_coefficients, fs::path(config$outputs$tables_dir, "abcd_sensitivity_coefficients.csv"))
safe_write_csv(abcd_results$sensitivity_full, fs::path(config$outputs$tables_dir, "abcd_sensitivity_full.csv"))
safe_write_csv(abcd_results$sensitivity_tests, fs::path(config$outputs$tables_dir, "abcd_sensitivity_tests.csv"))
safe_write_csv(abcd_results$precision_proxy_sensitivity, fs::path(config$outputs$tables_dir, "abcd_precision_proxy_sensitivity.csv"))

safe_write_csv(mcs_results$metrics, fs::path(config$outputs$tables_dir, "mcs_model_metrics.csv"))
safe_write_csv(mcs_results$coefficients, fs::path(config$outputs$tables_dir, "mcs_model_coefficients.csv"))
safe_write_csv(mcs_results$context_moderation, fs::path(config$outputs$tables_dir, "mcs_context_moderation.csv"))
safe_write_csv(mcs_results$context_full, fs::path(config$outputs$tables_dir, "mcs_context_full.csv"))
safe_write_csv(mcs_results$prs_moderation, fs::path(config$outputs$tables_dir, "mcs_prs_moderation.csv"))
safe_write_csv(mcs_results$prs_full, fs::path(config$outputs$tables_dir, "mcs_prs_full.csv"))
safe_write_csv(mcs_results$coefficient_tests, fs::path(config$outputs$tables_dir, "mcs_coefficient_tests.csv"))
safe_write_csv(mcs_results$sensitivity_coefficients, fs::path(config$outputs$tables_dir, "mcs_sensitivity_coefficients.csv"))
safe_write_csv(mcs_results$sensitivity_full, fs::path(config$outputs$tables_dir, "mcs_sensitivity_full.csv"))
safe_write_csv(mcs_results$sensitivity_tests, fs::path(config$outputs$tables_dir, "mcs_sensitivity_tests.csv"))
safe_write_csv(mcs_results$precision_proxy_sensitivity, fs::path(config$outputs$tables_dir, "mcs_precision_proxy_sensitivity.csv"))

### DT --> Missing-data sensitivity analyses for the primary models.
mcs_mi_primary <- fit_mcs_primary_mi_sensitivity(
  mcs,
  cache_path = fs::path(config$outputs$derived_dir, "mcs_primary_mi_m20_incomequintile_pttype2.rds"),
  m = 20
)
abcd_fiml_primary <- fit_abcd_primary_fiml_sensitivity(abcd)

safe_write_csv(mcs_mi_primary$coefficients, fs::path(config$outputs$tables_dir, "mcs_primary_mi_coefficients.csv"))
safe_write_csv(mcs_mi_primary$coefficient_test, fs::path(config$outputs$tables_dir, "mcs_primary_mi_test.csv"))
safe_write_csv(abcd_fiml_primary$coefficients, fs::path(config$outputs$tables_dir, "abcd_primary_fiml_coefficients.csv"))
safe_write_csv(abcd_fiml_primary$coefficient_test, fs::path(config$outputs$tables_dir, "abcd_primary_fiml_test.csv"))

### DT ---> MCS-parity sensitivity: ABCD primary model without parent education
abcd_no_parented <- fit_abcd_no_parented_sensitivity(abcd)
safe_write_csv(abcd_no_parented$coefficients,     fs::path(config$outputs$tables_dir, "abcd_no_parented_coefficients.csv"))
safe_write_csv(abcd_no_parented$coefficient_test, fs::path(config$outputs$tables_dir, "abcd_no_parented_test.csv"))

### DT --> ================================================================
### DT --> Balance-Index Analyses (Notebook 03)
### DT --> ================================================================
### DT --> Theory: ψ = log(distress / wellbeing-or-motivation) is an
### DT --> empirical proxy for the product γΔc — the combined location
### DT --> in the phase space governing transitions.
### DT --> MCS: ψ₁₇ = log(Kessler / WEMWBS), age 17 → age 23 outcomes
### DT --> ABCD: ψ = log(CBCL stress / Positive Affect), ses-03A → ses-05A
### DT --> ================================================================

source("R/balance_index.R")

mcs_balance_results <- fit_mcs_balance_models(mcs)
abcd_balance_results <- fit_abcd_balance_models(abcd, config)

safe_write_csv(mcs_balance_results$descriptives, fs::path(config$outputs$tables_dir, "mcs_balance_descriptives.csv"))
safe_write_csv(mcs_balance_results$comparison, fs::path(config$outputs$tables_dir, "mcs_balance_comparison.csv"))
safe_write_csv(mcs_balance_results$aic, fs::path(config$outputs$tables_dir, "mcs_balance_aic.csv"))

safe_write_csv(abcd_balance_results$descriptives, fs::path(config$outputs$tables_dir, "abcd_balance_descriptives.csv"))
safe_write_csv(abcd_balance_results$comparison, fs::path(config$outputs$tables_dir, "abcd_balance_comparison.csv"))
safe_write_csv(abcd_balance_results$aic, fs::path(config$outputs$tables_dir, "abcd_balance_aic.csv"))

message("Balance-index analyses complete.")

### DT --> ================================================================

claim_summary <- build_claim_summary(abcd_results, mcs_results)
safe_write_csv(claim_summary, fs::path(config$outputs$tables_dir, "claim_summary.csv"))

trajectory_abcd_plot <- make_trajectory_plot(
  abcd$long %>% filter(!is.na(internalising_t), session_id <= config$analysis$abcd_followup_session),
  "session_id",
  "internalising_t",
  "ABCD",
  "ABCD Internalising Trajectory Through Age 15"
)

trajectory_mcs_plot <- make_trajectory_plot(
  mcs$long %>% filter(!is.na(distress_z)),
  "wave",
  "distress_z",
  "MCS",
  "MCS Descriptive Distress Trajectory"
)

comparison_abcd_plot <- make_comparison_plot(
  abcd_results$coefficients,
  "ABCD",
  "ABCD Motivation/Engagement versus Cognitive-Control"
)

comparison_mcs_plot <- make_comparison_plot(
  mcs_results$coefficients,
  "MCS",
  "MCS Motivation/Engagement versus Executive-Control"
)

comparison_df <- bind_rows(
  comparison_plot_df(mcs_results$coefficients, "Millennium Cohort Study"),
  comparison_plot_df(abcd_results$coefficients, "Adolescent Brain Cognitive Development Study")
)
comparison_limits <- comparison_df %>%
  mutate(
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  ) %>%
  summarise(
    ymin = min(lower, 0, na.rm = TRUE) - 0.02,
    ymax = max(upper, 0, na.rm = TRUE) + 0.02
  )
comparison_ylim <- c(comparison_limits$ymin, comparison_limits$ymax)

comparison_abcd_plot <- make_comparison_plot(
  abcd_results$coefficients,
  "ABCD",
  "ABCD Motivation/Engagement versus Cognitive-Control",
  y_limits = comparison_ylim
)

comparison_mcs_plot <- make_comparison_plot(
  mcs_results$coefficients,
  "MCS",
  "MCS Motivation/Engagement versus Executive-Control",
  y_limits = comparison_ylim
)

comparison_combined_plot <- make_combined_comparison_plot(
  abcd_results,
  mcs_results,
  y_limits = comparison_ylim
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "abcd_trajectory.png"),
  plot = trajectory_abcd_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "mcs_trajectory.png"),
  plot = trajectory_mcs_plot,
  width = 8,
  height = 5,
  dpi = 300
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "abcd_field_vs_precision.png"),
  plot = comparison_abcd_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "abcd_field_vs_precision.pdf"),
  plot = comparison_abcd_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "mcs_field_vs_precision.png"),
  plot = comparison_mcs_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "mcs_field_vs_precision.pdf"),
  plot = comparison_mcs_plot,
  width = 7,
  height = 5
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "field_vs_precision_combined.png"),
  plot = comparison_combined_plot,
  width = 7.2,
  height = 3.6,
  dpi = 300
)

ggsave(
  filename = fs::path(config$outputs$figures_dir, "field_vs_precision_combined.pdf"),
  plot = comparison_combined_plot,
  width = 7.2,
  height = 3.6
)

summary_gt <- claim_summary %>% gt::gt()
safe_gtsave(summary_gt, fs::path(config$outputs$tables_dir, "claim_summary.html"))

message("Pipeline complete. Outputs written to ", fs::path(getwd(), "outputs"))
