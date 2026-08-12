### DT --> Re-run only the ABCD FIML sensitivity model from the existing
### DT --> analysis-ready cache. This script does not access data_ABCD.

source("R/utils.R")
source("R/models.R")

config <- load_project_config()
cache_path <- fs::path(config$outputs$derived_dir, "abcd_analysis_objects.rds")

if (!fs::file_exists(cache_path)) {
  stop("Analysis-ready ABCD cache not found: ", cache_path)
}

abcd <- readRDS(cache_path)
result <- fit_abcd_primary_fiml_sensitivity(abcd)

coefficient_path <- fs::path(
  config$outputs$tables_dir,
  "abcd_primary_fiml_coefficients.csv"
)
test_path <- fs::path(
  config$outputs$tables_dir,
  "abcd_primary_fiml_test.csv"
)

safe_write_csv(result$coefficients, coefficient_path)
safe_write_csv(result$coefficient_test, test_path)

contrast <- result$coefficient_test
if (
  nrow(contrast) != 1L ||
    !is.finite(contrast$std_error[[1]]) ||
    contrast$std_error[[1]] <= 0 ||
    !is.finite(contrast$statistic[[1]]) ||
    !is.finite(contrast$p_value[[1]])
) {
  stop("The regenerated FIML contrast failed numerical validation.")
}

reproduced_p <- 2 * stats::pnorm(
  abs(contrast$statistic[[1]]),
  lower.tail = FALSE
)
if (!isTRUE(all.equal(reproduced_p, contrast$p_value[[1]], tolerance = 1e-10))) {
  stop("The regenerated FIML z statistic and p value are inconsistent.")
}

message("ABCD FIML sensitivity re-run complete.")
message("Coefficient output: ", coefficient_path)
message("Contrast output: ", test_path)
message(
  sprintf(
    "Contrast estimate = %.6f; SE = %.6f; z = %.6f; p = %.6f; N = %d",
    contrast$estimate_difference[[1]],
    contrast$std_error[[1]],
    contrast$statistic[[1]],
    contrast$p_value[[1]],
    contrast$nobs[[1]]
  )
)
