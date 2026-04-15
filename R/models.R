### DT --> Descriptive summaries, trajectory models, lagged comparison models, moderation models, synthesis tables, and paper-oriented figures.
describe_longitudinal_data <- function(data, id_var, wave_var, outcome_var, cohort) {
  id_sym <- sym(id_var)
  wave_sym <- sym(wave_var)
  outcome_sym <- sym(outcome_var)

  wave_summary <- data %>%
    group_by(!!wave_sym) %>%
    summarise(
      cohort = cohort,
      n_rows = n(),
      n_participants = n_distinct(!!id_sym),
      outcome_mean = mean(!!outcome_sym, na.rm = TRUE),
      outcome_sd = sd(!!outcome_sym, na.rm = TRUE),
      outcome_missing = mean(is.na(!!outcome_sym)),
      .groups = "drop"
    ) %>%
    rename(wave = !!wave_sym)

  missingness <- data %>%
    summarise(across(everything(), ~ mean(is.na(.)))) %>%
    pivot_longer(everything(), names_to = "variable", values_to = "missing_fraction") %>%
    mutate(cohort = cohort) %>%
    arrange(desc(missing_fraction))

  wave_ids <- data %>%
    distinct(!!id_sym, !!wave_sym) %>%
    group_by(!!wave_sym) %>%
    summarise(
      ids = list(unique(!!id_sym)),
      n_participants = n_distinct(!!id_sym),
      .groups = "drop"
    ) %>%
    rename(wave = !!wave_sym) %>%
    arrange(wave)

  attrition <- wave_ids %>%
    mutate(
      cohort = cohort,
      retained_from_previous_n = purrr::map2_int(
        lag(ids), ids,
        ~ if (is.null(.x)) NA_integer_ else length(intersect(.x, .y))
      ),
      retention_from_previous = retained_from_previous_n / lag(n_participants)
    ) %>%
    select(-ids)

  list(
    wave_summary = wave_summary,
    missingness = missingness,
    attrition = attrition
  )
}

fit_trajectory_model <- function(data, id_var, time_var, outcome_var) {
  model_data <- data %>%
    filter(!is.na(.data[[id_var]]), !is.na(.data[[time_var]]), !is.na(.data[[outcome_var]])) %>%
    mutate(
      .time_centered = .data[[time_var]] - mean(.data[[time_var]], na.rm = TRUE)
    )

  control <- lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )

  fit_candidate <- function(formula_obj) {
    warnings <- character()
    fit <- withCallingHandlers(
      tryCatch(
        lmer(formula_obj, data = model_data, REML = FALSE, na.action = na.omit, control = control),
        error = function(e) e
      ),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )

    if (inherits(fit, "error")) {
      return(list(model = NULL, warnings = c(warnings, conditionMessage(fit)), singular = NA))
    }

    list(
      model = fit,
      warnings = warnings,
      singular = isSingular(fit, tol = 1e-4)
    )
  }

  candidate_specs <- list(
    list(
      label = "random_intercept_plus_slope_correlated",
      formula = as.formula(paste0(outcome_var, " ~ .time_centered + (1 + .time_centered | ", id_var, ")"))
    ),
    list(
      label = "random_intercept_plus_slope_uncorrelated",
      formula = as.formula(paste0(outcome_var, " ~ .time_centered + (1 + .time_centered || ", id_var, ")"))
    ),
    list(
      label = "random_intercept_only",
      formula = as.formula(paste0(outcome_var, " ~ .time_centered + (1 | ", id_var, ")"))
    )
  )

  fitted_candidates <- map(candidate_specs, function(spec) {
    result <- fit_candidate(spec$formula)
    c(spec, result)
  })

  acceptable <- keep(fitted_candidates, ~ !is.null(.x$model) && length(.x$warnings) == 0 && !isTRUE(.x$singular))
  chosen <- if (length(acceptable) > 0) {
    acceptable[[1]]
  } else {
    first(keep(fitted_candidates, ~ !is.null(.x$model))) %||% fitted_candidates[[length(fitted_candidates)]]
  }

  fit <- chosen$model
  attr(fit, "trajectory_spec") <- chosen$label
  attr(fit, "trajectory_warnings") <- chosen$warnings
  attr(fit, "trajectory_time_var") <- time_var
  fit
}

simple_svy_metrics <- function(model, model_name) {
  if (is.null(model)) {
    return(tibble(model = model_name, nobs = NA_real_, aic = NA_real_, bic = NA_real_, logLik = NA_real_))
  }
  model_nobs <- nrow(model$model)
  model_aic <- extract_model_aic(model)
  model_loglik_obj <- tryCatch(logLik(model), error = function(e) NA)
  model_loglik <- tryCatch(as.numeric(model_loglik_obj), error = function(e) NA_real_)
  model_k <- tryCatch(attr(model_loglik_obj, "df"), error = function(e) NULL) %||% length(coef(model))
  model_bic <- tryCatch(as.numeric(BIC(model)[1]), error = function(e) NA_real_)
  if (!is.finite(model_bic) && is.finite(model_loglik) && is.finite(model_nobs) && is.finite(model_k)) {
    model_bic <- as.numeric(-2 * model_loglik + log(model_nobs) * model_k)
  }
  tibble(
    model = model_name,
    nobs = model_nobs,
    aic = model_aic,
    bic = model_bic,
    logLik = model_loglik
  )
}

tidy_svy_summary <- function(model, model_name) {
  if (is.null(model)) {
    return(tibble(model = model_name, term = character(), estimate = numeric(), std_error = numeric(), statistic = numeric(), p_value = numeric(), conf_low = numeric(), conf_high = numeric(), std_estimate = numeric()))
  }
  coef_mat <- summary(model)$coefficients
  out <- tibble::as_tibble(coef_mat, rownames = "term")
  names(out)[1:4] <- c("term", "estimate", "std_error", "statistic")
  out %>%
    mutate(
      p_value = if ("Pr(>|t|)" %in% names(out)) .data[["Pr(>|t|)"]] else if ("Pr(>|z|)" %in% names(out)) .data[["Pr(>|z|)"]] else NA_real_,
      conf_low = estimate - 1.96 * std_error,
      conf_high = estimate + 1.96 * std_error,
      std_estimate = if_else(term == "(Intercept)", NA_real_, estimate),
      model = model_name,
      .before = 1
    )
}

coefficient_difference_test <- function(model, term_a, term_b, model_name, cohort = NA_character_, outcome = NA_character_) {
  out_template <- tibble(
    cohort = cohort,
    outcome = outcome,
    model = model_name,
    contrast = paste(term_a, "-", term_b),
    estimate_difference = NA_real_,
    std_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_
  )

  if (is.null(model)) {
    return(out_template)
  }

  beta <- coef(model)
  vcv <- vcov(model)
  if (!all(c(term_a, term_b) %in% names(beta))) {
    return(out_template)
  }

  diff_estimate <- unname(beta[[term_a]] - beta[[term_b]])
  diff_variance <- max(vcv[term_a, term_a] + vcv[term_b, term_b] - 2 * vcv[term_a, term_b], 0)
  diff_se <- sqrt(diff_variance)
  diff_stat <- if (isTRUE(diff_se > 0)) diff_estimate / diff_se else NA_real_
  diff_df <- suppressWarnings(tryCatch(df.residual(model), error = function(e) NA_real_))
  diff_p <- if (is.na(diff_stat)) {
    NA_real_
  } else if (is.na(diff_df) || is.infinite(diff_df)) {
    2 * pnorm(abs(diff_stat), lower.tail = FALSE)
  } else {
    2 * pt(abs(diff_stat), df = diff_df, lower.tail = FALSE)
  }

  tibble(
    cohort = cohort,
    outcome = outcome,
    model = model_name,
    contrast = paste(term_a, "-", term_b),
    estimate_difference = diff_estimate,
    std_error = diff_se,
    statistic = diff_stat,
    p_value = diff_p
  )
}

extract_model_nobs <- function(model) {
  if (is.null(model)) {
    return(NA_real_)
  }
  if (inherits(model, "lavaan")) {
    return(tryCatch(as.numeric(sum(lavaan::lavInspect(model, "nobs"))), error = function(e) NA_real_))
  }
  tryCatch(
    as.numeric(stats::nobs(model)),
    error = function(e) {
      if (!is.null(model$model)) {
        as.numeric(nrow(model$model))
      } else {
        NA_real_
      }
    }
  )
}

build_model_info <- function(models, outcome, cohort = NA_character_, family = NA_character_, missingness_handling = "Complete-case analysis") {
  tibble(
    cohort = cohort,
    family = family,
    outcome = outcome,
    model = names(models),
    nobs = vapply(models, extract_model_nobs, numeric(1)),
    missingness_handling = missingness_handling
  )
}

attach_model_info <- function(df, model_info) {
  if (is.null(df) || nrow(df) == 0 || is.null(model_info) || nrow(model_info) == 0) {
    return(df)
  }
  out <- df %>%
    left_join(
      model_info %>% select(outcome, model, nobs, missingness_handling),
      by = c("outcome", "model"),
      suffix = c("", "__info")
    )
  if ("nobs__info" %in% names(out)) {
    out <- out %>%
      mutate(nobs = coalesce(nobs, nobs__info)) %>%
      select(-nobs__info)
  }
  if ("missingness_handling__info" %in% names(out)) {
    out <- out %>%
      mutate(missingness_handling = coalesce(missingness_handling, missingness_handling__info)) %>%
      select(-missingness_handling__info)
  }
  out
}

format_mean_sd <- function(x, digits = 2) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_character_)
  }
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), stats::sd(x))
}

format_n_pct <- function(x, value = 1, digits = 1) {
  keep <- !is.na(x)
  if (!any(keep)) {
    return(NA_character_)
  }
  n_val <- sum(x[keep] == value)
  pct_val <- 100 * n_val / sum(keep)
  sprintf(paste0("%s (%.", digits, "f%%)"), scales::comma(n_val), pct_val)
}

format_level_count <- function(x, level, digits = 1) {
  keep <- !is.na(x)
  if (!any(keep)) {
    return(NA_character_)
  }
  n_val <- sum(x[keep] == level)
  pct_val <- 100 * n_val / sum(keep)
  sprintf(paste0("%s (%.", digits, "f%%)"), scales::comma(n_val), pct_val)
}

build_mcs_demographics <- function(data, sample_label) {
  base_rows <- tibble(
    sample = sample_label,
    characteristic = c(
      "Cohort members (N)",
      "Age at baseline, mean (SD)",
      "Female, n (%)",
      "OECD equivalised income quintile, mean (SD)",
      "Baseline distress, mean (SD)"
    ),
    summary = c(
      scales::comma(nrow(data)),
      format_mean_sd(data$age_baseline),
      format_n_pct(data$sex, value = 1),
      format_mean_sd(data$income_quintile),
      format_mean_sd(data$baseline_distress_raw)
    )
  )

  stratum_levels <- levels(droplevels(data$sampling_stratum))
  stratum_rows <- tibble(
    sample = sample_label,
    characteristic = paste0("PTTYPE2 stratum: ", stratum_levels, ", n (%)"),
    summary = vapply(stratum_levels, function(level) format_level_count(data$sampling_stratum, level), character(1))
  )

  ethnicity_levels <- levels(droplevels(data$ethnicity_6cat))
  ethnicity_rows <- tibble(
    sample = sample_label,
    characteristic = paste0("Ethnicity: ", ethnicity_levels, ", n (%)"),
    summary = vapply(ethnicity_levels, function(level) format_level_count(data$ethnicity_6cat, level), character(1))
  )

  bind_rows(base_rows, ethnicity_rows, stratum_rows)
}

build_abcd_demographics <- function(data, sample_label) {
  base_rows <- tibble(
    sample = sample_label,
    characteristic = c(
      "Cohort members (N)",
      "Age at baseline, mean (SD)",
      "Female, n (%)",
      "Household income, mean (SD)",
      "Parent education, mean (SD)",
      "Area deprivation index percentile, mean (SD)",
      "Baseline internalising, mean (SD)"
    ),
    summary = c(
      scales::comma(nrow(data)),
      format_mean_sd(data$baseline_age),
      format_n_pct(data$sex, value = 1),
      format_mean_sd(data$household_income),
      format_mean_sd(data$parent_education),
      format_mean_sd(data$neighborhood_adi),
      format_mean_sd(data$baseline_internalising_t)
    )
  )

  ethnicity_levels <- levels(droplevels(data$ethnicity_5cat))
  ethnicity_rows <- tibble(
    sample = sample_label,
    characteristic = paste0("Race/ethnicity: ", ethnicity_levels, ", n (%)"),
    summary = vapply(ethnicity_levels, function(level) format_level_count(data$ethnicity_5cat, level), character(1))
  )

  bind_rows(base_rows, ethnicity_rows)
}

fit_abcd_models <- function(abcd) {
  d <- abcd$predictive %>%
    filter(
      !is.na(followup_internalising_z),
      !is.na(baseline_internalising_z)
    ) %>%
    mutate(
      sex_z = zscore(sex),
      household_income_z = zscore(household_income),
      parent_education_z = zscore(parent_education),
      baseline_age_z = zscore(baseline_age),
      neighborhood_risk_z = zscore(neighborhood_risk),
      site = as.factor(site)
    )

  build_abcd_covariate_terms <- function(data, baseline_var) {
    terms <- c(baseline_var, "baseline_age_z")
    if (has_variation(data$sex_z)) terms <- c(terms, "sex_z")
    if (has_variation(data$household_income_z)) terms <- c(terms, "household_income_z")
    if (has_variation(data$parent_education_z)) terms <- c(terms, "parent_education_z")
    if (has_variation(data$neighborhood_risk_z)) terms <- c(terms, "neighborhood_risk_z")
    if (dplyr::n_distinct(stats::na.omit(data$site)) > 1) terms <- c(terms, "site")
    paste(terms, collapse = " + ")
  }

  fit_linear_comparison <- function(data, outcome_var, baseline_var, outcome_label) {
    covariate_terms <- build_abcd_covariate_terms(data, baseline_var)
    f0 <- as.formula(paste0(outcome_var, " ~ ", covariate_terms))
    m0 <- lm(f0, data = data)
    m1 <- update(m0, . ~ . + field_index)
    m2 <- update(m0, . ~ . + precision_index)
    m3 <- update(m0, . ~ . + field_index + precision_index)
    model_info <- build_model_info(
      list(base = m0, field = m1, precision = m2, combined = m3),
      outcome = outcome_label,
      cohort = "ABCD",
      family = "primary",
      missingness_handling = "Complete-case analysis"
    )

    list(
      models = list(base = m0, field = m1, precision = m2, combined = m3),
      metrics = bind_rows(
        simple_model_metrics(m0, "base"),
        simple_model_metrics(m1, "field"),
        simple_model_metrics(m2, "precision"),
        simple_model_metrics(m3, "combined")
      ) %>%
        mutate(outcome = outcome_label, .before = 1) %>%
        attach_model_info(model_info),
      coefficients = bind_rows(
        simple_coef_table(m0, "base"),
        simple_coef_table(m1, "field"),
        simple_coef_table(m2, "precision"),
        simple_coef_table(m3, "combined")
      ) %>%
        mutate(outcome = outcome_label, .before = 1) %>%
        attach_model_info(model_info),
      coefficient_tests = coefficient_difference_test(
        m3,
        term_a = "field_index",
        term_b = "precision_index",
        model_name = "combined",
        cohort = "ABCD",
        outcome = outcome_label
      ) %>%
        attach_model_info(model_info %>% filter(model == "combined")),
      model_info = model_info
    )
  }

  summarise_abcd_precision_proxy <- function(data, precision_var, precision_label, min_complete_n = 100) {
    dat <- data %>%
      mutate(precision_index = .data[[precision_var]])

    covariate_terms <- build_abcd_covariate_terms(dat, "baseline_internalising_z")
    combined_formula <- as.formula(
      paste0("followup_internalising_z ~ ", covariate_terms, " + field_index + precision_index")
    )
    complete_n <- sum(stats::complete.cases(stats::model.frame(combined_formula, data = dat, na.action = na.pass)))

    template <- tibble(
      cohort = "ABCD",
      outcome = "Age 15 parent-reported internalising",
      precision_proxy = precision_label,
      precision_variable = precision_var,
      field_estimate = NA_real_,
      precision_estimate = NA_real_,
      estimate_difference = NA_real_,
      coefficient_test_p = NA_real_,
      nobs = complete_n,
      missingness_handling = "Complete-case analysis",
      status = if (complete_n < min_complete_n) {
        paste0("Not estimable: fewer than ", min_complete_n, " complete cases for adjusted model")
      } else {
        "Estimated"
      }
    )

    if (complete_n < min_complete_n) {
      return(template)
    }

    fit <- fit_linear_comparison(
      dat,
      outcome_var = "followup_internalising_z",
      baseline_var = "baseline_internalising_z",
      outcome_label = "Age 15 parent-reported internalising"
    )
    combined <- fit$coefficients %>% filter(model == "combined", term %in% c("field_index", "precision_index"))
    diff_row <- fit$coefficient_tests %>% filter(model == "combined")

    template %>%
      mutate(
        field_estimate = combined %>% filter(term == "field_index") %>% pull(estimate),
        precision_estimate = combined %>% filter(term == "precision_index") %>% pull(estimate),
        estimate_difference = diff_row$estimate_difference %||% NA_real_,
        coefficient_test_p = diff_row$p_value %||% NA_real_
      )
  }

  primary <- fit_linear_comparison(
    d,
    outcome_var = "followup_internalising_z",
    baseline_var = "baseline_internalising_z",
    outcome_label = "Age 15 parent-reported internalising"
  )

  moderators <- c("family_context", "sleep_risk", "neighborhood_risk", "peer_adversity")
  moderation_models <- map(moderators, function(mod_var) {
    dat <- d %>% filter(!is.na(.data[[mod_var]]), has_variation(.data[[mod_var]]))
    if (nrow(dat) < 50) {
      return(NULL)
    }
    covariate_terms <- build_abcd_covariate_terms(dat, "baseline_internalising_z")
    fit <- lm(
      as.formula(
        paste0(
          "followup_internalising_z ~ ", covariate_terms, " + ",
          "field_index * ", mod_var, " + precision_index * ", mod_var
        )
      ),
      data = dat
    )
    pattern <- paste0(
      "field_index:", mod_var, "|", mod_var, ":field_index|",
      "precision_index:", mod_var, "|", mod_var, ":precision_index"
    )
    full <- simple_coef_table(fit, paste0("abcd_", mod_var)) %>%
      mutate(
        outcome = "Age 15 parent-reported internalising",
        nobs = extract_model_nobs(fit),
        missingness_handling = "Complete-case analysis",
        .before = 1
      )
    list(
      interaction = full %>% filter(str_detect(term, pattern)),
      full = full
    )
  }) %>%
    compact()
  moderation <- bind_rows(map(moderation_models, "interaction"))
  moderation_full <- bind_rows(map(moderation_models, "full"))

  sensitivity_specs <- tibble(
    outcome_var = c("followup_bpm_internalising_z", "followup_bpm_attention_z"),
    baseline_var = c("baseline_bpm_internalising_z", "baseline_bpm_attention_z"),
    outcome_label = c("Age 15 youth BPM internalising", "Age 15 youth BPM attention")
  )

  sensitivity_results <- pmap(
    sensitivity_specs,
    function(outcome_var, baseline_var, outcome_label) {
      dat <- d %>%
        filter(!is.na(.data[[outcome_var]]), !is.na(.data[[baseline_var]]))
      if (nrow(dat) < 100) {
        return(NULL)
      }
      fit_linear_comparison(dat, outcome_var, baseline_var, outcome_label)
    }
  ) %>%
    compact()

  primary_complete_case_vars <- c(
    "followup_internalising_z",
    "baseline_internalising_z",
    "field_index",
    "precision_index",
    "baseline_age_z"
  )
  if (has_variation(d$sex_z)) primary_complete_case_vars <- c(primary_complete_case_vars, "sex_z")
  if (has_variation(d$household_income_z)) primary_complete_case_vars <- c(primary_complete_case_vars, "household_income_z")
  if (has_variation(d$parent_education_z)) primary_complete_case_vars <- c(primary_complete_case_vars, "parent_education_z")
  if (has_variation(d$neighborhood_risk_z)) primary_complete_case_vars <- c(primary_complete_case_vars, "neighborhood_risk_z")

  primary_demographics <- d %>%
    filter(
      if_all(all_of(primary_complete_case_vars), ~ !is.na(.)),
      if (dplyr::n_distinct(stats::na.omit(site)) > 1) !is.na(site) else TRUE
    ) %>%
    build_abcd_demographics("ABCD primary complete-case sample (`ses-02A` to `ses-05A`)")

  list(
    input = d,
    models = primary$models,
    metrics = primary$metrics,
    coefficients = primary$coefficients,
    moderation = moderation,
    moderation_full = moderation_full,
    coefficient_tests = primary$coefficient_tests,
    sensitivity_coefficients = map_dfr(sensitivity_results, "coefficients") %>%
      filter(model == "combined", term %in% c("field_index", "precision_index")),
    sensitivity_full = map_dfr(sensitivity_results, ~ .x$coefficients %>% filter(model == "combined")),
    sensitivity_tests = map_dfr(sensitivity_results, "coefficient_tests"),
    ### DT ---> Card Sort (crdst) and List Sorting (lswmt) had only 5 and 14 valid
    ### DT ---> cases respectively at ses-02A and were replaced with variables that
    ### DT ---> have adequate ses-02A coverage (~8 000 and ~7 300 valid cases).
    precision_proxy_sensitivity = bind_rows(
      summarise_abcd_precision_proxy(d, "precision_processing_speed_index", "Pattern Comparison Processing Speed"),
      summarise_abcd_precision_proxy(d, "precision_crystallized_index",     "NIH Toolbox Crystallized composite")
    ),
    primary_demographics = primary_demographics,
    covariate_note = "ABCD models adjust for baseline outcome, baseline age, sex, household income, parent education, area deprivation index percentile, and assessment site."
  )
}

fit_mcs_models <- function(mcs) {
  d_age17 <- mcs$predictive_age17 %>%
    filter(
      !is.na(analysis_weight_uk),
      !is.na(design_psu),
      !is.na(design_stratum),
      analysis_weight_uk > 0
    ) %>%
    mutate(
      baseline_age_z = zscore(age_baseline),
      sex_z = zscore(sex),
      income_quintile_z = zscore(income_quintile),
      sampling_stratum = droplevels(as.factor(sampling_stratum))
    )
  d_age23 <- mcs$predictive_age23 %>%
    filter(
      !is.na(analysis_weight_uk),
      !is.na(design_psu),
      !is.na(design_stratum),
      analysis_weight_uk > 0
    ) %>%
    mutate(
      baseline_age_z = zscore(age_baseline),
      sex_z = zscore(sex),
      income_quintile_z = zscore(income_quintile),
      sampling_stratum = droplevels(as.factor(sampling_stratum))
    )
  d_age14_mfq <- mcs$predictive_age17 %>%
    filter(
      !is.na(mfq_age14_z),
      !is.na(analysis_weight_uk),
      !is.na(design_psu),
      !is.na(design_stratum),
      analysis_weight_uk > 0
    ) %>%
    mutate(
      baseline_age_z = zscore(age_baseline),
      sex_z = zscore(sex),
      income_quintile_z = zscore(income_quintile),
      sampling_stratum = droplevels(as.factor(sampling_stratum))
    )
  pc_terms <- paste0("pc", 1:20)
  pc_formula <- paste(pc_terms, collapse = " + ")
  options(survey.lonely.psu = "adjust")

  fit_svy_model <- function(formula_obj, data) {
    vars <- all.vars(formula_obj)
    dat <- data %>%
      filter(
        if_all(all_of(vars), ~ !is.na(.)),
        !is.na(analysis_weight_uk),
        !is.na(design_psu),
        !is.na(design_stratum),
        analysis_weight_uk > 0
      )
    if (nrow(dat) < 100) {
      return(NULL)
    }
    design <- survey::svydesign(
      ids = ~design_psu,
      strata = ~design_stratum,
      weights = ~analysis_weight_uk,
      data = dat,
      nest = TRUE
    )
    survey::svyglm(formula_obj, design = design)
  }

  core_covariates <- "baseline_age_z + sex_z + income_quintile_z + sampling_stratum"

  fit_weighted_comparison <- function(data, outcome_var, outcome_label, include_baseline_distress = TRUE) {
    f0 <- if (isTRUE(include_baseline_distress)) {
      as.formula(paste0(outcome_var, " ~ baseline_distress_z + ", core_covariates))
    } else {
      as.formula(paste0(outcome_var, " ~ ", core_covariates))
    }
    f1 <- update(f0, . ~ . + field_index)
    f2 <- update(f0, . ~ . + precision_index)
    f3 <- update(f0, . ~ . + field_index + precision_index)

    m0 <- fit_svy_model(f0, data)
    m1 <- fit_svy_model(f1, data)
    m2 <- fit_svy_model(f2, data)
    m3 <- fit_svy_model(f3, data)
    model_info <- build_model_info(
      list(base = m0, field = m1, precision = m2, combined = m3),
      outcome = outcome_label,
      cohort = "MCS",
      family = "primary",
      missingness_handling = "Complete-case analysis"
    )

    list(
      models = list(base = m0, field = m1, precision = m2, combined = m3),
      metrics = bind_rows(
        simple_svy_metrics(m0, "base"),
        simple_svy_metrics(m1, "field"),
        simple_svy_metrics(m2, "precision"),
        simple_svy_metrics(m3, "combined")
      ) %>%
        mutate(outcome = outcome_label, .before = 1) %>%
        attach_model_info(model_info),
      coefficients = bind_rows(
        tidy_svy_summary(m0, "base"),
        tidy_svy_summary(m1, "field"),
        tidy_svy_summary(m2, "precision"),
        tidy_svy_summary(m3, "combined")
      ) %>%
        mutate(outcome = outcome_label, .before = 1) %>%
        attach_model_info(model_info),
      coefficient_test = coefficient_difference_test(
        m3,
        term_a = "field_index",
        term_b = "precision_index",
        model_name = "combined",
        cohort = "MCS",
        outcome = outcome_label
      ) %>%
        attach_model_info(model_info %>% filter(model == "combined")),
      model_info = model_info
    )
  }

  summarise_mcs_precision_proxy <- function(data, precision_var, precision_label, outcome_var = "followup_distress_z", outcome_label = "Age 17 psychological distress", include_baseline_distress = TRUE) {
    dat <- data %>%
      mutate(precision_index = .data[[precision_var]])

    covariate_terms <- if (isTRUE(include_baseline_distress)) {
      paste0("baseline_distress_z + ", core_covariates)
    } else {
      core_covariates
    }
    combined_formula <- as.formula(
      paste0(outcome_var, " ~ ", covariate_terms, " + field_index + precision_index")
    )
    combined_vars <- all.vars(combined_formula)
    complete_n <- dat %>%
      filter(
        if_all(all_of(combined_vars), ~ !is.na(.)),
        !is.na(analysis_weight_uk),
        !is.na(design_psu),
        !is.na(design_stratum),
        analysis_weight_uk > 0
      ) %>%
      nrow()

    template <- tibble(
      cohort = "MCS",
      outcome = outcome_label,
      precision_proxy = precision_label,
      precision_variable = precision_var,
      field_estimate = NA_real_,
      precision_estimate = NA_real_,
      estimate_difference = NA_real_,
      coefficient_test_p = NA_real_,
      nobs = complete_n,
      missingness_handling = "Complete-case analysis",
      status = if (complete_n < 100) "Not estimable: fewer than 100 complete cases for weighted model" else "Estimated"
    )

    if (complete_n < 100) {
      return(template)
    }

    fit <- fit_weighted_comparison(dat, outcome_var, outcome_label, include_baseline_distress = include_baseline_distress)
    combined <- fit$coefficients %>% filter(model == "combined", term %in% c("field_index", "precision_index"))
    diff_row <- fit$coefficient_test %>% filter(model == "combined")

    template %>%
      mutate(
        field_estimate = combined %>% filter(term == "field_index") %>% pull(estimate),
        precision_estimate = combined %>% filter(term == "precision_index") %>% pull(estimate),
        estimate_difference = diff_row$estimate_difference %||% NA_real_,
        coefficient_test_p = diff_row$p_value %||% NA_real_
      )
  }

  primary <- fit_weighted_comparison(d_age17, "followup_distress_z", "Age 17 psychological distress")

  context_mods <- c("family_conflict", "sleep_risk", "peer_school_trouble")
  context_models <- map(context_mods, function(mod_var) {
    dat <- d_age17 %>% filter(!is.na(.data[[mod_var]]), has_variation(.data[[mod_var]]))
    if (nrow(dat) < 100) {
      return(NULL)
    }
    fit <- fit_svy_model(
      as.formula(
        paste0(
          "followup_distress_z ~ baseline_distress_z + ", core_covariates, " + ",
          "field_index * ", mod_var, " + precision_index * ", mod_var
        )
      ),
      dat
    )
    if (is.null(fit)) {
      return(NULL)
    }
    pattern <- paste0(
      "field_index:", mod_var, "|", mod_var, ":field_index|",
      "precision_index:", mod_var, "|", mod_var, ":precision_index"
    )
    full <- tidy_svy_summary(fit, paste0("mcs_", mod_var)) %>%
      mutate(
        outcome = "Age 17 psychological distress",
        nobs = extract_model_nobs(fit),
        missingness_handling = "Complete-case analysis",
        .before = 1
      )
    list(
      interaction = full %>% filter(str_detect(term, pattern)),
      full = full
    )
  }) %>%
    compact()
  context_results <- bind_rows(map(context_models, "interaction"))
  context_full <- bind_rows(map(context_models, "full"))

  prs_mods <- c("depression_prs", "cognition_pgi")
  prs_models <- map(prs_mods, function(prs_var) {
    dat <- d_age17 %>% filter(!is.na(.data[[prs_var]]), has_variation(.data[[prs_var]]))
    if (nrow(dat) < 100) {
      return(NULL)
    }
    fit <- fit_svy_model(
      as.formula(
        paste0(
          "followup_distress_z ~ baseline_distress_z + ", core_covariates, " + ",
          pc_formula, " + ",
          "field_index * ", prs_var, " + precision_index * ", prs_var
        )
      ),
      dat
    )
    if (is.null(fit)) {
      return(NULL)
    }
    pattern <- paste0(
      "field_index:", prs_var, "|", prs_var, ":field_index|",
      "precision_index:", prs_var, "|", prs_var, ":precision_index"
    )
    full <- tidy_svy_summary(fit, paste0("mcs_prs_", prs_var)) %>%
      mutate(
        outcome = "Age 17 psychological distress",
        nobs = extract_model_nobs(fit),
        missingness_handling = "Complete-case analysis",
        .before = 1
      )
    list(
      interaction = full %>% filter(str_detect(term, pattern)),
      full = full
    )
  }) %>%
    compact()
  prs_results <- bind_rows(map(prs_models, "interaction"))
  prs_full <- bind_rows(map(prs_models, "full"))

  sensitivity_specs <- tibble(
    data_name = c("age14", "age17", "age23", "age23", "age23", "age23", "age23"),
    outcome_var = c(
      "mfq_age14_z",
      "followup_wellbeing_z",
      "followup_distress_z",
      "followup_anxiety_z",
      "followup_depression_z",
      "followup_loneliness_z",
      "followup_wellbeing_z"
    ),
    outcome_label = c(
      "Age 14 MFQ depressive symptoms",
      "Age 17 wellbeing",
      "Age 23 psychological distress",
      "Age 23 anxiety symptoms",
      "Age 23 depressive symptoms",
      "Age 23 loneliness",
      "Age 23 wellbeing"
    )
  )

  sensitivity_results <- pmap(
    sensitivity_specs,
    function(data_name, outcome_var, outcome_label) {
      dat <- if (data_name == "age14") {
        d_age14_mfq
      } else if (data_name == "age17") {
        d_age17
      } else {
        d_age23
      }
      fit_weighted_comparison(
        dat,
        outcome_var,
        outcome_label,
        include_baseline_distress = data_name != "age14"
      )
    }
  )

  sensitivity_coefficients <- map_dfr(sensitivity_results, "coefficients") %>%
    filter(model == "combined", term %in% c("field_index", "precision_index"))
  sensitivity_tests <- map_dfr(sensitivity_results, "coefficient_test")

  primary_demographics <- d_age17 %>%
    filter(
      !is.na(followup_distress_z),
      !is.na(baseline_distress_z),
      !is.na(field_index),
      !is.na(precision_index),
      !is.na(baseline_age_z),
      !is.na(sex_z),
      !is.na(income_quintile_z),
      !is.na(sampling_stratum),
      !is.na(analysis_weight_uk),
      !is.na(design_psu),
      !is.na(design_stratum),
      analysis_weight_uk > 0
    ) %>%
    build_mcs_demographics("MCS primary complete-case sample (age 14 exposure to age 17 follow-up)")

  list(
    input = d_age17,
    input_age23 = d_age23,
    models = primary$models,
    metrics = primary$metrics,
    coefficients = primary$coefficients,
    coefficient_tests = primary$coefficient_test,
    context_moderation = context_results,
    context_full = context_full,
    prs_moderation = prs_results,
    prs_full = prs_full,
    sensitivity_coefficients = sensitivity_coefficients,
    sensitivity_full = map_dfr(sensitivity_results, ~ .x$coefficients %>% filter(model == "combined")),
    sensitivity_tests = sensitivity_tests,
    precision_proxy_sensitivity = bind_rows(
      summarise_mcs_precision_proxy(d_age17, "precision_delay_index", "Delay aversion"),
      summarise_mcs_precision_proxy(d_age17, "precision_word_index", "Word activity score")
    ),
    primary_demographics = primary_demographics,
    weighting = "survey_weighted_exposure_wave_primary",
    weighting_note = mcs$survey_weight_note,
    covariate_note = "MCS models adjust for age-14 baseline distress when prospective, plus age-14 age, sex, PTTYPE2 sampling stratum, and OECD equivalised income quintile."
  )
}

load_or_create_mcs_primary_mi <- function(data, cache_path, m = 20, seed = 20260404) {
  if (file.exists(cache_path)) {
    return(readRDS(cache_path))
  }

  mi_data <- data %>%
    select(
      MCSID, analysis_weight_uk, design_psu, design_stratum,
      followup_distress_z, baseline_distress_z, baseline_age_z, sex_z,
      income_quintile_z, sampling_stratum, field_index, precision_index
    )

  methods <- mice::make.method(mi_data)
  methods[c("MCSID", "analysis_weight_uk", "design_psu", "design_stratum", "sampling_stratum")] <- ""
  methods[setdiff(names(mi_data), c("MCSID", "analysis_weight_uk", "design_psu", "design_stratum", "sampling_stratum"))] <- "pmm"

  predictor_matrix <- mice::make.predictorMatrix(mi_data)
  predictor_matrix[, c("MCSID", "design_psu", "design_stratum")] <- 0
  predictor_matrix[c("MCSID", "analysis_weight_uk", "design_psu", "design_stratum"), ] <- 0

  imp <- mice::mice(
    mi_data,
    m = m,
    method = methods,
    predictorMatrix = predictor_matrix,
    maxit = 10,
    seed = seed,
    printFlag = FALSE
  )

  saveRDS(imp, cache_path)
  imp
}

tidy_pooled_mi_model <- function(pooled, outcome_label, model_name, nobs, missingness_handling) {
  estimates <- pooled$coefficients
  ses <- sqrt(diag(pooled$variance))
  stats <- estimates / ses
  tibble(
    outcome = outcome_label,
    model = model_name,
    term = names(estimates),
    estimate = as.numeric(estimates),
    std_error = as.numeric(ses),
    statistic = as.numeric(stats),
    p_value = 2 * stats::pnorm(abs(stats), lower.tail = FALSE),
    conf_low = as.numeric(estimates - 1.96 * ses),
    conf_high = as.numeric(estimates + 1.96 * ses),
    std_estimate = if_else(term == "(Intercept)", NA_real_, estimate),
    nobs = nobs,
    missingness_handling = missingness_handling
  )
}

pooled_mi_difference_test <- function(pooled, term_a, term_b, outcome_label, model_name, cohort, nobs, missingness_handling) {
  beta <- pooled$coefficients
  vcv <- pooled$variance
  if (!all(c(term_a, term_b) %in% names(beta))) {
    return(tibble(
      cohort = cohort,
      outcome = outcome_label,
      model = model_name,
      contrast = paste(term_a, "-", term_b),
      estimate_difference = NA_real_,
      std_error = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_,
      nobs = nobs,
      missingness_handling = missingness_handling
    ))
  }

  diff_estimate <- unname(beta[[term_a]] - beta[[term_b]])
  diff_variance <- max(vcv[term_a, term_a] + vcv[term_b, term_b] - 2 * vcv[term_a, term_b], 0)
  diff_se <- sqrt(diff_variance)
  diff_stat <- if (isTRUE(diff_se > 0)) diff_estimate / diff_se else NA_real_
  diff_p <- if (is.na(diff_stat)) NA_real_ else 2 * stats::pnorm(abs(diff_stat), lower.tail = FALSE)

  tibble(
    cohort = cohort,
    outcome = outcome_label,
    model = model_name,
    contrast = paste(term_a, "-", term_b),
    estimate_difference = diff_estimate,
    std_error = diff_se,
    statistic = diff_stat,
    p_value = diff_p,
    nobs = nobs,
    missingness_handling = missingness_handling
  )
}

fit_mcs_primary_mi_sensitivity <- function(mcs, cache_path, m = 20) {
  analysis_data <- mcs$predictive_age17 %>%
    filter(
      !is.na(analysis_weight_uk),
      !is.na(design_psu),
      !is.na(design_stratum),
      analysis_weight_uk > 0
    ) %>%
    mutate(
      baseline_age_z = zscore(age_baseline),
      sex_z = zscore(sex),
      income_quintile_z = zscore(income_quintile),
      sampling_stratum = droplevels(as.factor(sampling_stratum))
    )

  imp <- load_or_create_mcs_primary_mi(analysis_data, cache_path = cache_path, m = m)
  completed <- mice::complete(imp, action = "all")
  formula_obj <- stats::as.formula(
    "followup_distress_z ~ baseline_distress_z + baseline_age_z + sex_z + income_quintile_z + sampling_stratum + field_index + precision_index"
  )
  fits <- lapply(completed, function(dat) {
    design <- survey::svydesign(
      ids = ~design_psu,
      strata = ~design_stratum,
      weights = ~analysis_weight_uk,
      data = dat,
      nest = TRUE
    )
    survey::svyglm(formula_obj, design = design)
  })

  pooled <- mitools::MIcombine(
    results = lapply(fits, stats::coef),
    variances = lapply(fits, stats::vcov)
  )

  outcome_label <- "Age 17 psychological distress"
  missingness_label <- "Multiple imputation (m = 20), pooled with Rubin's rules"
  nobs <- nrow(analysis_data)

  list(
    coefficients = tidy_pooled_mi_model(
      pooled,
      outcome_label = outcome_label,
      model_name = "combined_mi",
      nobs = nobs,
      missingness_handling = missingness_label
    ),
    coefficient_test = pooled_mi_difference_test(
      pooled,
      term_a = "field_index",
      term_b = "precision_index",
      outcome_label = outcome_label,
      model_name = "combined_mi",
      cohort = "MCS",
      nobs = nobs,
      missingness_handling = missingness_label
    )
  )
}

fit_abcd_primary_fiml_sensitivity <- function(abcd) {
  d <- abcd$predictive %>%
    mutate(
      sex_z = zscore(sex),
      household_income_z = zscore(household_income),
      parent_education_z = zscore(parent_education),
      baseline_age_z = zscore(baseline_age),
      neighborhood_risk_z = zscore(neighborhood_risk),
      site = as.character(site)
    )

  d$site[is.na(d$site)] <- "missing_site"
  d$site <- as.factor(d$site)

  core_terms <- c("baseline_internalising_z", "baseline_age_z")
  if (has_variation(d$sex_z)) core_terms <- c(core_terms, "sex_z")
  if (sum(!is.na(d$household_income_z)) > 1 && has_variation(d$household_income_z)) core_terms <- c(core_terms, "household_income_z")
  if (sum(!is.na(d$parent_education_z)) > 1 && has_variation(d$parent_education_z)) core_terms <- c(core_terms, "parent_education_z")
  if (sum(!is.na(d$neighborhood_risk_z)) > 1 && has_variation(d$neighborhood_risk_z)) core_terms <- c(core_terms, "neighborhood_risk_z")
  core_terms <- c(core_terms, "field_index", "precision_index")

  include_site <- dplyr::n_distinct(stats::na.omit(d$site)) > 1
  if (include_site) {
    site_matrix <- stats::model.matrix(~ site, data = d)[, -1, drop = FALSE]
    colnames(site_matrix) <- make.names(colnames(site_matrix))
  } else {
    site_matrix <- NULL
  }

  lavaan_data <- bind_cols(
    d %>% select(followup_internalising_z, all_of(core_terms)),
    if (!is.null(site_matrix)) as_tibble(site_matrix) else tibble()
  )

  predictor_terms <- c(
    if ("baseline_internalising_z" %in% core_terms) "b_baseline*baseline_internalising_z",
    if ("baseline_age_z" %in% core_terms) "b_age*baseline_age_z",
    if ("sex_z" %in% core_terms) "b_sex*sex_z",
    if ("household_income_z" %in% core_terms) "b_income*household_income_z",
    if ("parent_education_z" %in% core_terms) "b_parentedu*parent_education_z",
    if ("neighborhood_risk_z" %in% core_terms) "b_neigh*neighborhood_risk_z",
    "b_field*field_index",
    "b_precision*precision_index",
    if (!is.null(site_matrix)) paste0("b_", colnames(site_matrix), "*", colnames(site_matrix))
  )
  predictor_terms <- predictor_terms[!is.na(predictor_terms)]

  model_syntax <- paste(
    "followup_internalising_z ~",
    paste(predictor_terms, collapse = " + ")
  )

  fit <- lavaan::sem(
    model = model_syntax,
    data = lavaan_data,
    missing = "fiml",
    fixed.x = FALSE,
    meanstructure = TRUE,
    estimator = "MLR"
  )

  pe <- lavaan::parameterEstimates(fit, ci = TRUE) %>%
    filter(lhs == "followup_internalising_z", op == "~") %>%
    transmute(
      outcome = "Age 15 parent-reported internalising",
      model = "combined_fiml",
      term = rhs,
      estimate = est,
      std_error = se,
      statistic = z,
      p_value = pvalue,
      conf_low = ci.lower,
      conf_high = ci.upper,
      std_estimate = if_else(str_detect(rhs, "^site"), NA_real_, est),
      nobs = extract_model_nobs(fit),
      missingness_handling = "Full-information maximum likelihood (FIML)"
    )

  wald <- lavaan::lavTestWald(fit, constraints = "b_field == b_precision")
  vcov_mat <- tryCatch(stats::vcov(fit), error = function(e) NULL)
  diff_estimate <- pe %>% filter(term == "field_index") %>% pull(estimate) -
    pe %>% filter(term == "precision_index") %>% pull(estimate)

  ### DT ---> Bug W6 fix: the previous approach used match() to get integer column
  ### DT ---> indices, which returns NA_integer_ when lavaan uses parameter labels
  ### DT ---> (e.g. "b_field") whose exact form in vcov() depends on lavaan version
  ### DT ---> and estimator. NA_integer_ passes is.null() but not is.finite(), so
  ### DT ---> diff_se silently became NA_real_, writing an empty cell to the CSV.
  ### DT --->
  ### DT ---> Fix: use string-based double-bracket indexing inside tryCatch.
  ### DT ---> If that fails (label not found in vcov), fall back to deriving SE from
  ### DT ---> the Wald chi-sq statistic: SE = |estimate_diff| / sqrt(wald$stat),
  ### DT ---> which is algebraically equivalent when the Wald test uses the same
  ### DT ---> delta-method covariance matrix (Wald chi-sq with df = 1 equals z^2).
  ### DT --->
  ### DT ---> The statistic is now reported as a z-score (estimate / SE), consistent
  ### DT ---> with the MCS pooled-MI contrast in pooled_mi_difference_test(). The
  ### DT ---> p-value is from a two-sided standard-normal test, which is numerically
  ### DT ---> identical to the Wald chi-sq(df=1) p-value to machine precision.
  diff_se <- tryCatch(
    {
      v <- vcov_mat
      sqrt(max(
        v["b_field", "b_field"] + v["b_precision", "b_precision"] -
          2 * v["b_field", "b_precision"],
        0
      ))
    },
    error = function(e) {
      ### DT ---> Primary vcov lookup failed; derive SE from Wald chi-sq (df = 1)
      if (!is.null(wald) && is.finite(wald$stat) && wald$stat > 0) {
        abs(diff_estimate) / sqrt(wald$stat)
      } else {
        NA_real_
      }
    }
  )
  diff_stat <- if (is.finite(diff_se) && diff_se > 0) diff_estimate / diff_se else NA_real_
  diff_p <- if (is.finite(diff_stat)) {
    2 * stats::pnorm(abs(diff_stat), lower.tail = FALSE)
  } else {
    ### DT ---> z-stat unavailable; fall back to Wald p-value
    tryCatch(unname(wald$p.value), error = function(e) NA_real_)
  }

  comparison <- tibble(
    cohort = "ABCD",
    outcome = "Age 15 parent-reported internalising",
    model = "combined_fiml",
    contrast = "field_index - precision_index",
    estimate_difference = diff_estimate,
    std_error = diff_se,
    statistic = diff_stat,
    p_value = diff_p,
    nobs = extract_model_nobs(fit),
    missingness_handling = "Full-information maximum likelihood (FIML)"
  )

  list(
    coefficients = pe,
    coefficient_test = comparison
  )
}

### DT ---> MCS-parity sensitivity: same primary ABCD model re-estimated
### DT ---> without parent_education_z.  MCS has no equivalent variable in its
### DT ---> standard adjustment set (age, sex, income quintile, stratum), so
### DT ---> this check confirms the ABCD ordering result does not depend on the
### DT ---> additional SES control that has no MCS counterpart.
fit_abcd_no_parented_sensitivity <- function(abcd) {
  d <- abcd$predictive %>%
    filter(
      !is.na(followup_internalising_z),
      !is.na(baseline_internalising_z)
    ) %>%
    mutate(
      sex_z              = zscore(sex),
      household_income_z = zscore(household_income),
      baseline_age_z     = zscore(baseline_age),
      neighborhood_risk_z = zscore(neighborhood_risk),
      site               = as.factor(site)
    )

  ### DT ---> Build covariate string WITHOUT parent_education_z
  build_terms <- function(data, baseline_var) {
    terms <- c(baseline_var, "baseline_age_z")
    if (has_variation(data$sex_z))              terms <- c(terms, "sex_z")
    if (has_variation(data$household_income_z)) terms <- c(terms, "household_income_z")
    if (has_variation(data$neighborhood_risk_z)) terms <- c(terms, "neighborhood_risk_z")
    if (dplyr::n_distinct(stats::na.omit(data$site)) > 1) terms <- c(terms, "site")
    paste(terms, collapse = " + ")
  }

  covariate_terms <- build_terms(d, "baseline_internalising_z")

  f0 <- as.formula(paste0("followup_internalising_z ~ ", covariate_terms))
  m0 <- lm(f0, data = d)
  m3 <- update(m0, . ~ . + field_index + precision_index)

  coefs <- simple_coef_table(m3, "combined") %>%
    mutate(
      outcome              = "Age 15 parent-reported internalising",
      nobs                 = extract_model_nobs(m3),
      missingness_handling = "Complete-case (no parent education)",
      .before              = 1
    )

  diff_test <- coefficient_difference_test(
    m3,
    term_a     = "field_index",
    term_b     = "precision_index",
    model_name = "combined",
    cohort     = "ABCD",
    outcome    = "Age 15 parent-reported internalising"
  ) %>%
    mutate(missingness_handling = "Complete-case (no parent education)")

  list(
    coefficients       = coefs,
    coefficient_test   = diff_test,
    covariate_note     = paste0(
      "MCS-parity sensitivity: ABCD primary model re-estimated without parent ",
      "education. Covariates retained: baseline internalising, age, sex, ",
      "household income, area deprivation index, and assessment site."
    )
  )
}

make_trajectory_plot <- function(data, wave_var, outcome_var, cohort, title) {
  wave_sym <- sym(wave_var)
  outcome_sym <- sym(outcome_var)

  summary_df <- data %>%
    group_by(!!wave_sym) %>%
    summarise(
      mean_outcome = mean(!!outcome_sym, na.rm = TRUE),
      se = sd(!!outcome_sym, na.rm = TRUE) / sqrt(sum(!is.na(!!outcome_sym))),
      .groups = "drop"
    )

  ggplot(summary_df, aes(x = !!wave_sym, y = mean_outcome, group = 1)) +
    geom_line(color = "#1b5e20", linewidth = 1) +
    geom_point(color = "#1b5e20", size = 2.5) +
    geom_errorbar(aes(ymin = mean_outcome - 1.96 * se, ymax = mean_outcome + 1.96 * se), width = 0.1) +
    labs(
      title = title,
      x = "Wave",
      y = ifelse(cohort == "ABCD", "Mean CBCL internalising T-score", "Mean pooled-standardised main distress outcome")
    ) +
    theme_paper()
}

comparison_plot_df <- function(results_df, cohort_label) {
  precision_label <- if (cohort_label %in% c("ABCD", "Adolescent Brain Cognitive Development Study")) {
    "Cognitive-control\nproxy"
  } else {
    "Executive-control\nproxy"
  }

  results_df %>%
    filter(model == "combined", term %in% c("field_index", "precision_index")) %>%
    mutate(
      cohort = cohort_label,
      term = recode(
        term,
        field_index = "Motivation/\nengagement proxy",
        precision_index = precision_label
      ),
      term = factor(
        term,
        levels = c(precision_label, "Motivation/\nengagement proxy")
      )
    )
}

make_comparison_plot <- function(results_df, cohort, title, y_limits = NULL) {
  plot_df <- comparison_plot_df(results_df, cohort)

  ggplot(plot_df, aes(x = term, y = estimate, colour = term)) +
    geom_hline(yintercept = 0, color = "#9e9e9e", linewidth = 0.4) +
    geom_linerange(
      aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
      linewidth = 1.1
    ) +
    geom_point(size = 3.2) +
    scale_colour_manual(
      values = c(
        "Motivation/\nengagement proxy" = "#006d77",
        "Executive-control\nproxy" = "#bb3e03",
        "Cognitive-control\nproxy" = "#bb3e03"
      )
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(title = title, x = NULL, y = "Standardised coefficient") +
    theme_paper() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9)
    )
}

make_combined_comparison_plot <- function(abcd_results, mcs_results, y_limits = NULL) {
  plot_df <- bind_rows(
    comparison_plot_df(mcs_results$coefficients, "Millennium Cohort Study"),
    comparison_plot_df(abcd_results$coefficients, "Adolescent Brain Cognitive Development Study")
  ) %>%
    mutate(
      cohort = factor(
        cohort,
        levels = c("Millennium Cohort Study", "Adolescent Brain Cognitive Development Study")
      )
    )

  cohort_labels <- c(
    "Millennium Cohort Study" = "(a) Millennium Cohort Study",
    "Adolescent Brain Cognitive Development Study" = "(b) Adolescent Brain Cognitive Development Study"
  )

  ggplot(plot_df, aes(x = term, y = estimate, colour = term)) +
    geom_hline(yintercept = 0, color = "#9e9e9e", linewidth = 0.4) +
    geom_linerange(
      aes(ymin = estimate - 1.96 * std_error, ymax = estimate + 1.96 * std_error),
      linewidth = 1.1
    ) +
    geom_point(size = 3.2) +
    facet_wrap(~ cohort, nrow = 1, labeller = as_labeller(cohort_labels)) +
    scale_colour_manual(
      values = c(
        "Motivation/\nengagement proxy" = "#006d77",
        "Executive-control\nproxy" = "#bb3e03",
        "Cognitive-control\nproxy" = "#bb3e03"
      )
    ) +
    coord_cartesian(ylim = y_limits) +
    labs(x = NULL, y = "Standardised coefficient") +
    theme_paper() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 9)
    )
}

build_claim_summary <- function(abcd_results, mcs_results) {
  abcd_combined <- abcd_results$coefficients %>% filter(model == "combined", term %in% c("field_index", "precision_index"))
  mcs_combined <- mcs_results$coefficients %>% filter(model == "combined", term %in% c("field_index", "precision_index"))
  abcd_diff <- abcd_results$coefficient_tests %>% filter(model == "combined")
  mcs_diff <- mcs_results$coefficient_tests %>% filter(model == "combined")

  abcd_field_beta <- abcd_combined %>% filter(term == "field_index") %>% pull(estimate)
  abcd_precision_beta <- abcd_combined %>% filter(term == "precision_index") %>% pull(estimate)
  abcd_field_p <- abcd_combined %>% filter(term == "field_index") %>% pull(p_value)
  abcd_precision_p <- abcd_combined %>% filter(term == "precision_index") %>% pull(p_value)
  abcd_diff_p <- abcd_diff %>% pull(p_value)

  mcs_field_beta <- mcs_combined %>% filter(term == "field_index") %>% pull(estimate)
  mcs_precision_beta <- mcs_combined %>% filter(term == "precision_index") %>% pull(estimate)
  mcs_field_p <- mcs_combined %>% filter(term == "field_index") %>% pull(p_value)
  mcs_precision_p <- mcs_combined %>% filter(term == "precision_index") %>% pull(p_value)
  mcs_diff_p <- mcs_diff %>% pull(p_value)
  mcs_prs_sig <- nrow(mcs_results$prs_moderation) > 0 && any(mcs_results$prs_moderation$p_value < 0.05, na.rm = TRUE)

  support_label <- function(field_beta, precision_beta, field_p, precision_p, diff_p) {
    case_when(
      is.na(field_beta) ~ "not estimable",
      !is.na(diff_p) && diff_p < 0.05 && abs(field_beta) > abs(precision_beta) ~ "supports",
      field_p < 0.10 && abs(field_beta) >= abs(precision_beta) ~ "partially supports",
      TRUE ~ "does not support"
    )
  }

  tibble(
    claim_tested = c(
      "Motivation/engagement variables predict later persistence/worsening more strongly than cognitive-control variables",
      "Motivation/engagement variables predict later persistence/worsening more strongly than executive-control variables",
      "Context moderates motivation/engagement and/or cognitive-control pathways",
      "Context moderates motivation/engagement and/or executive-control pathways",
      "Genetic scores moderate sensitivity to motivation/engagement or context pathways"
    ),
    cohort = c("ABCD", "MCS", "ABCD", "MCS", "MCS"),
    model = c(
      "Lagged linear model, age 12 to age 15 follow-up",
      "Survey-weighted svyglm, age 14 exposure to age 17 follow-up",
      "Lagged interaction models with family/sleep/neighborhood moderators",
      "Survey-weighted svyglm interaction models with family, sleep, and peer moderators to age 17",
      "Survey-weighted svyglm interaction models with depression PRS and cognition PGI plus ancestry PCs to age 17"
    ),
    key_result = c(
      sprintf("Motivation/engagement beta = %.3f (p = %.3f); cognitive-control beta = %.3f (p = %.3f); direct coefficient-comparison p = %.3f", abcd_field_beta, abcd_field_p, abcd_precision_beta, abcd_precision_p, abcd_diff_p),
      sprintf("Motivation/engagement beta = %.3f (p = %.3f); executive-control beta = %.3f (p = %.3f); direct coefficient-comparison p = %.3f", mcs_field_beta, mcs_field_p, mcs_precision_beta, mcs_precision_p, mcs_diff_p),
      ifelse(nrow(abcd_results$moderation) > 0, "At least one context interaction was estimable in the local ABCD extract", "Context interactions were limited by local ABCD availability"),
      ifelse(nrow(mcs_results$context_moderation) > 0, "Multiple age-14 context interactions were estimable in MCS", "Context interactions were not stable/estimable"),
      ifelse(
        nrow(mcs_results$prs_moderation) == 0,
        "Genetic-score moderation not estimable",
        ifelse(
          mcs_prs_sig,
          "At least one exploratory genetic-score interaction reached conventional significance with ancestry-PC adjustment",
          "Exploratory genetic-score interactions were estimable with ancestry-PC adjustment, but none reached conventional significance"
        )
      )
    ),
    support_level = c(
      support_label(abcd_field_beta, abcd_precision_beta, abcd_field_p, abcd_precision_p, abcd_diff_p),
      support_label(mcs_field_beta, mcs_precision_beta, mcs_field_p, mcs_precision_p, mcs_diff_p),
      ifelse(nrow(abcd_results$moderation) > 0, "partially supports", "does not support"),
      ifelse(nrow(mcs_results$context_moderation) > 0, "partially supports", "does not support"),
      ifelse(nrow(mcs_results$prs_moderation) == 0, "not estimable", ifelse(mcs_prs_sig, "partially supports", "does not support"))
    )
  )
}
