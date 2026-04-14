### DT --> Shared helpers for loading local cohort data, standardising scales, and writing outputs.
suppressPackageStartupMessages({
  library(dplyr)
  library(fs)
  library(ggplot2)
  library(gt)
  library(haven)
  library(jsonlite)
  library(lme4)
  library(lmerTest)
  library(purrr)
  library(readr)
  library(readxl)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) {
    y
  } else {
    x
  }
}

load_project_config <- function(path = "config/project_config.yml") {
  yaml::read_yaml(path)
}

project_path <- function(...) {
  fs::path(getwd(), ...)
}

ensure_directories <- function(paths) {
  walk(paths, fs::dir_create, recurse = TRUE)
}

clean_numeric <- function(x) {
  if (inherits(x, "haven_labelled")) {
    x <- haven::zap_labels(x)
  }
  if (inherits(x, "POSIXt")) {
    return(x)
  }
  if (is.character(x)) {
    x <- na_if(x, "n/a")
    x <- na_if(x, "")
    suppressWarnings(x <- as.numeric(x))
  }
  x <- as.numeric(x)
  x[x %in% c(444, 555, 666, 777, 888, 999)] <- NA_real_
  x
}

clean_character <- function(x) {
  x %>%
    na_if("n/a") %>%
    na_if("")
}

zscore <- function(x) {
  x <- clean_numeric(x)
  if (sum(!is.na(x)) < 2 || isTRUE(sd(x, na.rm = TRUE) == 0)) {
    return(rep(NA_real_, length(x)))
  }
  as.numeric(scale(x))
}

row_mean_standardised <- function(data, vars) {
  if (length(vars) == 0) {
    return(rep(NA_real_, nrow(data)))
  }
  z_mat <- data %>%
    mutate(across(all_of(vars), zscore)) %>%
    select(all_of(vars)) %>%
    as.matrix()
  rowMeans(z_mat, na.rm = TRUE) %>%
    replace(is.nan(.), NA_real_)
}

row_sum_with_min <- function(data, vars, min_non_missing = length(vars), transform = identity) {
  if (length(vars) == 0) {
    return(rep(NA_real_, nrow(data)))
  }
  mat <- data %>%
    transmute(across(all_of(vars), ~ transform(clean_numeric(.)))) %>%
    as.matrix()
  non_missing <- rowSums(!is.na(mat))
  out <- rowSums(mat, na.rm = TRUE)
  out[non_missing < min_non_missing] <- NA_real_
  out
}

first_non_missing <- function(...) {
  vectors <- list(...)
  reduce(vectors, coalesce)
}

first_non_missing_scalar <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA)
  }
  x[[1]]
}

collapse_by_id <- function(data, id_var) {
  id_sym <- sym(id_var)
  data %>%
    group_by(!!id_sym) %>%
    summarise(across(everything(), first_non_missing_scalar), .groups = "drop")
}

read_abcd_tsv <- function(path, cols = NULL) {
  readr::read_tsv(
    file = path,
    col_select = all_of(cols %||% character()),
    na = c("", "n/a", "NA"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(across(where(is.character), clean_character))
}

read_abcd_tsv_any <- function(path, cols) {
  readr::read_tsv(
    file = path,
    col_select = any_of(cols),
    na = c("", "n/a", "NA"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(across(where(is.character), clean_character))
}

read_abcd_parquet_any <- function(path, cols) {
  if (!fs::is_absolute_path(path)) {
    path <- fs::path(getwd(), path)
  }
  if (!file.exists(path)) {
    stop("ABCD parquet file not found: ", path)
  }

  tmp_out <- tempfile(fileext = ".tsv")
  tmp_py <- tempfile(fileext = ".py")
  tmp_log <- tempfile(fileext = ".log")
  py_code <- c(
    "import sys",
    "import pyarrow.parquet as pq",
    "import pyarrow.csv as pcsv",
    "",
    "path = sys.argv[1]",
    "requested = [c for c in sys.argv[2].split('__COLSEP__') if c]",
    "schema = pq.read_schema(path).names",
    "cols = [c for c in requested if c in schema]",
    "table = pq.read_table(path, columns=cols)",
    "pcsv.write_csv(table, sys.argv[3], write_options=pcsv.WriteOptions(delimiter='\\t'))"
  )
  writeLines(py_code, tmp_py)

  python_candidates <- unique(c(
    "/opt/anaconda3/bin/python3",
    Sys.which("python3"),
    "/usr/bin/python3"
  ))
  python_candidates <- python_candidates[nzchar(python_candidates) & file.exists(python_candidates)]

  usable_python <- NULL
  if (length(python_candidates) > 0) {
    for (candidate in python_candidates) {
      probe_status <- suppressWarnings(system2(
        candidate,
        args = c("-c", "import pyarrow.parquet as pq"),
        stdout = FALSE,
        stderr = FALSE
      ))
      if (identical(as.integer(probe_status), 0L)) {
        usable_python <- candidate
        break
      }
    }
  }

  if (is.null(usable_python) || !nzchar(usable_python)) {
    stop("Could not locate a usable python3 binary for parquet import.")
  }

  status <- system2(
    usable_python,
    args = c(tmp_py, path, paste(cols, collapse = "__COLSEP__"), tmp_out),
    stdout = tmp_log,
    stderr = tmp_log
  )

  if (!identical(as.integer(status), 0L) || !file.exists(tmp_out)) {
    log_text <- if (file.exists(tmp_log)) paste(readLines(tmp_log, warn = FALSE), collapse = "\n") else ""
    stop(
      "Failed to read ABCD parquet file: ", path,
      if (nzchar(log_text)) paste0("\nPython log:\n", log_text) else ""
    )
  }

  on.exit(unlink(c(tmp_out, tmp_py, tmp_log)), add = TRUE)

  readr::read_tsv(
    file = tmp_out,
    col_select = any_of(cols),
    na = c("", "n/a", "NA"),
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    mutate(across(where(is.character), clean_character))
}

read_mcs_sav_any <- function(path, cols) {
  haven::read_sav(file = path, col_select = any_of(cols))
}

filter_mcs_member_number <- function(data, member_var, target_value) {
  if (!member_var %in% names(data)) {
    stop("Required MCS member identifier not found: ", member_var)
  }
  data %>%
    filter(clean_numeric(.data[[member_var]]) == target_value)
}

filter_mcs_first_cm <- function(data, member_var) {
  filter_mcs_member_number(data, member_var, target_value = 1)
}

filter_mcs_first_pnum <- function(data, member_var = "pnum") {
  filter_mcs_member_number(data, member_var, target_value = 100)
}

session_to_wave_year <- function(session_id) {
  case_when(
    session_id == "ses-00A" ~ 0,
    session_id == "ses-01A" ~ 1,
    session_id == "ses-02A" ~ 2,
    session_id == "ses-03A" ~ 3,
    session_id == "ses-04A" ~ 4,
    session_id == "ses-05A" ~ 5,
    session_id == "ses-06A" ~ 6,
    TRUE ~ NA_real_
  )
}

theme_paper <- function() {
  theme_minimal(base_size = 10, base_family = "sans") +
    theme(
      plot.background = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      strip.background = element_rect(fill = "white", colour = NA),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "#e6e6e6", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.35),
      axis.ticks = element_line(color = "black", linewidth = 0.35),
      axis.ticks.length = grid::unit(2, "pt"),
      plot.title = element_text(face = "bold", size = 10, hjust = 0),
      axis.title = element_text(face = "bold"),
      strip.text = element_text(face = "bold", size = 10)
    )
}

safe_write_csv <- function(data, path) {
  readr::write_csv(data, path, na = "")
  invisible(path)
}

safe_write_rds <- function(object, path) {
  saveRDS(object, path)
  invisible(path)
}

safe_gtsave <- function(data, path) {
  tryCatch(
    {
      gt::gtsave(data, path)
      path
    },
    error = function(e) {
      warning("Could not save gt output to ", path, ": ", conditionMessage(e))
      NULL
    }
  )
}

extract_model_aic <- function(model) {
  if (is.null(model)) {
    return(NA_real_)
  }
  raw_aic <- tryCatch(extractAIC(model), error = function(e) NA_real_)
  if (length(raw_aic) == 0 || all(is.na(raw_aic))) {
    return(NA_real_)
  }
  if (!is.null(names(raw_aic)) && "AIC" %in% names(raw_aic)) {
    return(as.numeric(raw_aic[["AIC"]]))
  }
  if (length(raw_aic) >= 2) {
    return(as.numeric(raw_aic[[2]]))
  }
  as.numeric(raw_aic[[1]])
}

simple_model_metrics <- function(model, model_name) {
  fitted_model <- model
  model_nobs <- tryCatch(
    stats::nobs(fitted_model),
    error = function(e) length(stats::residuals(fitted_model))
  )
  tibble(
    model = model_name,
    nobs = model_nobs,
    aic = extract_model_aic(fitted_model),
    bic = BIC(fitted_model),
    logLik = as.numeric(logLik(fitted_model))
  )
}

simple_coef_table <- function(model, model_name) {
  coef_mat <- summary(model)$coefficients
  tibble::as_tibble(coef_mat, rownames = "term") %>%
    rename_with(~ c("estimate", "std_error", "statistic", "p_value")[seq_along(.)], -term) %>%
    mutate(
      conf_low = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error,
      std_estimate = if_else(term == "(Intercept)", NA_real_, estimate),
      model = model_name,
      .before = 1
    )
}

tidy_lmer_summary <- function(model, model_name) {
  coef_mat <- summary(model)$coefficients
  tibble::as_tibble(coef_mat, rownames = "term") %>%
    mutate(term = if_else(term == ".time_centered", attr(model, "trajectory_time_var") %||% term, term)) %>%
    rename(
      estimate = Estimate,
      std_error = `Std. Error`,
      statistic = `t value`
    ) %>%
    mutate(
      conf_low = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error,
      std_estimate = if_else(term == "(Intercept)", NA_real_, estimate),
      p_value = suppressWarnings(2 * pt(abs(statistic), df = summary(model)$devcomp$dims[["n"]] - 1, lower.tail = FALSE)),
      trajectory_spec = attr(model, "trajectory_spec") %||% NA_character_,
      model = model_name,
      .before = 1
    )
}

has_variation <- function(x) {
  x <- x[!is.na(x)]
  length(unique(x)) > 1
}
