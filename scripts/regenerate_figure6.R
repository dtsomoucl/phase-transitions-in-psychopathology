### DT --> Regenerate the main cohort comparison figure from non-disclosive
### DT --> aggregate coefficient tables. This script does not access raw data.

source("R/utils.R")
source("R/models.R")

config <- load_project_config()
table_dir <- config$outputs$tables_dir
figure_dir <- config$outputs$figures_dir

abcd_coefficients <- readr::read_csv(
  fs::path(table_dir, "abcd_model_coefficients.csv"),
  show_col_types = FALSE
)
mcs_coefficients <- readr::read_csv(
  fs::path(table_dir, "mcs_model_coefficients.csv"),
  show_col_types = FALSE
)

combined_coefficients <- dplyr::bind_rows(
  comparison_plot_df(mcs_coefficients, "Millennium Cohort Study"),
  comparison_plot_df(
    abcd_coefficients,
    "Adolescent Brain Cognitive Development Study"
  )
)

comparison_limits <- combined_coefficients %>%
  dplyr::mutate(
    lower = estimate - 1.96 * std_error,
    upper = estimate + 1.96 * std_error
  ) %>%
  dplyr::summarise(
    ymin = min(lower, 0, na.rm = TRUE) - 0.02,
    ymax = max(upper, 0, na.rm = TRUE) + 0.02
  )

comparison_plot <- make_combined_comparison_plot(
  abcd_results = list(coefficients = abcd_coefficients),
  mcs_results = list(coefficients = mcs_coefficients),
  y_limits = c(comparison_limits$ymin, comparison_limits$ymax)
)

pdf_path <- fs::path(figure_dir, "field_vs_precision_combined.pdf")
png_path <- fs::path(figure_dir, "field_vs_precision_combined.png")
workspace_pdf_path <- fs::path("..", "Fig_6.pdf")
workspace_png_path <- fs::path("..", "Fig_6.png")

ggplot2::ggsave(pdf_path, comparison_plot, width = 7.2, height = 3.6)
ggplot2::ggsave(png_path, comparison_plot, width = 7.2, height = 3.6, dpi = 300)
file.copy(pdf_path, workspace_pdf_path, overwrite = TRUE)
file.copy(png_path, workspace_png_path, overwrite = TRUE)

message("Figure 6 regenerated from aggregate coefficient tables.")
message("Workspace PDF: ", workspace_pdf_path)
message("Workspace PNG: ", workspace_png_path)
