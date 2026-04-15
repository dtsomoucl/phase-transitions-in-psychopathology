### DT --> Balance-index analyses (Notebook 03 companion code)
### DT --> psi = log(distress / wellbeing-or-motivation) is treated here as a
### DT --> theory-derived observable intended to index proximity to a clinical
### DT --> transition. The outcome models below therefore use binary probable-
### DT --> case outcomes rather than continuous symptom scores.

balance_note_result <- function(msg) {
  list(
    descriptives = tibble::tibble(note = msg),
    comparison = tibble::tibble(note = msg),
    aic = tibble::tibble(note = msg)
  )
}

balance_aic_table <- function(comparison) {
  comparison %>%
    dplyr::select(outcome, model, aic, auc, brier, nobs, missingness_handling) %>%
    dplyr::distinct() %>%
    dplyr::group_by(outcome) %>%
    dplyr::mutate(
      delta_aic = if (any(is.finite(aic))) {
        aic - min(aic[is.finite(aic)], na.rm = TRUE)
      } else {
        rep(NA_real_, dplyr::n())
      }
    ) %>%
    dplyr::ungroup()
}

balance_model_aic <- function(model) {
  if (is.null(model)) {
    return(NA_real_)
  }
  if (inherits(model, "svyglm")) {
    raw_aic <- tryCatch(extractAIC(model), error = function(e) NA_real_)
    if (length(raw_aic) > 0 && any(is.finite(raw_aic))) {
      if (!is.null(names(raw_aic)) && "AIC" %in% names(raw_aic)) {
        return(as.numeric(raw_aic[["AIC"]]))
      }
      if (length(raw_aic) >= 2) {
        return(as.numeric(raw_aic[[2]]))
      }
      return(as.numeric(raw_aic[[1]]))
    }
    return(tryCatch(as.numeric(stats::deviance(model) + 2 * length(stats::coef(model))), error = function(e) NA_real_))
  }
  ll <- tryCatch(stats::logLik(model), error = function(e) NA)
  if (all(is.finite(ll))) {
    df_ll <- attr(ll, "df") %||% length(stats::coef(model))
    return(as.numeric(-2 * ll + 2 * df_ll))
  }
  tryCatch(as.numeric(stats::deviance(model) + 2 * length(stats::coef(model))), error = function(e) NA_real_)
}

balance_weighted_auc <- function(y, score, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, length(y))
  }

  df <- tibble::tibble(
    y = as.numeric(y),
    score = as.numeric(score),
    w = as.numeric(w)
  ) %>%
    dplyr::filter(
      is.finite(y),
      is.finite(score),
      is.finite(w),
      w > 0,
      y %in% c(0, 1)
    )

  if (nrow(df) == 0) {
    return(NA_real_)
  }

  total_pos <- sum(df$w[df$y == 1], na.rm = TRUE)
  total_neg <- sum(df$w[df$y == 0], na.rm = TRUE)
  if (!is.finite(total_pos) || !is.finite(total_neg) || total_pos <= 0 || total_neg <= 0) {
    return(NA_real_)
  }

  grouped <- df %>%
    dplyr::arrange(score) %>%
    dplyr::group_by(score) %>%
    dplyr::summarise(
      w_pos = sum(w[y == 1], na.rm = TRUE),
      w_neg = sum(w[y == 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(cum_neg_before = dplyr::lag(cumsum(w_neg), default = 0))

  auc_num <- sum(grouped$w_pos * (grouped$cum_neg_before + 0.5 * grouped$w_neg), na.rm = TRUE)
  auc_num / (total_pos * total_neg)
}

balance_brier_score <- function(y, prob, w = NULL) {
  if (is.null(w)) {
    w <- rep(1, length(y))
  }

  keep <- is.finite(y) & is.finite(prob) & is.finite(w) & w > 0
  if (!any(keep)) {
    return(NA_real_)
  }

  y <- as.numeric(y[keep])
  prob <- as.numeric(prob[keep])
  w <- as.numeric(w[keep])
  sum(w * (y - prob)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

balance_model_fit_stats <- function(model) {
  if (is.null(model)) {
    return(list(aic = NA_real_, auc = NA_real_, brier = NA_real_, nobs = NA_real_))
  }

  obs <- attr(model, "balance_obs")
  weights <- attr(model, "balance_weights")
  prob <- tryCatch(stats::predict(model, type = "response"), error = function(e) rep(NA_real_, length(obs)))
  nobs <- tryCatch(stats::nobs(model), error = function(e) length(obs))

  list(
    aic = balance_model_aic(model),
    auc = balance_weighted_auc(obs, prob, weights),
    brier = balance_brier_score(obs, prob, weights),
    nobs = nobs
  )
}

balance_extract_coef <- function(model, predictor, model_label, outcome_label) {
  if (is.null(model)) {
    return(NULL)
  }

  coef_mat <- summary(model)$coefficients
  if (!predictor %in% rownames(coef_mat)) {
    return(NULL)
  }

  stat_col <- intersect(c("t value", "z value"), colnames(coef_mat))[1]
  p_col <- intersect(c("Pr(>|t|)", "Pr(>|z|)"), colnames(coef_mat))[1]
  estimate <- unname(coef_mat[predictor, "Estimate"])
  std_error <- unname(coef_mat[predictor, "Std. Error"])
  fit_stats <- balance_model_fit_stats(model)

  tibble::tibble(
    outcome = outcome_label,
    model = model_label,
    predictor = predictor,
    estimate = estimate,
    std_error = std_error,
    statistic = if (!is.na(stat_col)) unname(coef_mat[predictor, stat_col]) else NA_real_,
    p_value = if (!is.na(p_col)) unname(coef_mat[predictor, p_col]) else NA_real_,
    odds_ratio = exp(estimate),
    or_conf_low = exp(estimate - 1.96 * std_error),
    or_conf_high = exp(estimate + 1.96 * std_error),
    aic = fit_stats$aic,
    auc = fit_stats$auc,
    brier = fit_stats$brier,
    nobs = fit_stats$nobs
  )
}

make_mcs_balance_subset <- function(df) {
  kessler_candidates <- c(
    "gdckessl_17", "GDCKESSL", "gdckessl", "kessler_17",
    "distress_17", "kessler6_17", "gdckessl_17_raw"
  )
  wemwbs_candidates <- c(
    "gdwemwbs_17", "GDWEMWBS", "gdwemwbs", "wemwbs_17",
    "wellbeing_17", "gdwemwbs_17_raw"
  )

  kessler_var <- intersect(kessler_candidates, names(df))
  wemwbs_var <- intersect(wemwbs_candidates, names(df))

  if (length(kessler_var) == 0 || length(wemwbs_var) == 0) {
    return(NULL)
  }

  kessler_var <- kessler_var[[1]]
  wemwbs_var <- wemwbs_var[[1]]

  bal <- df %>%
    dplyr::rename(
      kessler_val = !!rlang::sym(kessler_var),
      wemwbs_val = !!rlang::sym(wemwbs_var)
    ) %>%
    dplyr::filter(
      !is.na(kessler_val), !is.na(wemwbs_val),
      kessler_val > 0, wemwbs_val > 0
    ) %>%
    dplyr::mutate(
      age17_age = first_non_missing(age17_age, age_baseline + 3),
      age17_age_z = zscore(age17_age),
      psi_17 = log(kessler_val / wemwbs_val),
      psi_17_z = zscore(psi_17),
      kessler_17_z = zscore(kessler_val),
      wemwbs_17_z = zscore(wemwbs_val),
      depression_case_23 = dplyr::if_else(!is.na(depression_23), as.numeric(depression_23 >= 3), NA_real_),
      GOVWT2 = clean_numeric(GOVWT2),
      PTTYPE2 = clean_numeric(PTTYPE2),
      SPTN00 = clean_character(SPTN00)
    )

  if (!"sex_z" %in% names(bal) || all(is.na(bal$sex_z))) bal$sex_z <- zscore(bal$sex)
  if (!"income_quintile_z" %in% names(bal) || all(is.na(bal$income_quintile_z))) {
    income_source <- if ("income_quintile" %in% names(bal)) bal$income_quintile else bal$income_quintile_z
    bal$income_quintile_z <- zscore(income_source)
  }
  if (!"sampling_stratum" %in% names(bal) || all(is.na(bal$sampling_stratum))) bal$sampling_stratum <- mcs_pttype2_factor(bal$PTTYPE2)
  bal$sampling_stratum <- droplevels(as.factor(bal$sampling_stratum))

  bal
}

fit_mcs_balance_core <- function(bal, suffix = "") {
  if (is.null(bal) || nrow(bal) == 0) {
    return(list(
      descriptives = tibble::tibble(note = "No estimable models"),
      comparison = tibble::tibble(note = "No estimable models"),
      aic = tibble::tibble(note = "No estimable models")
    ))
  }

  desc <- bal %>%
    dplyr::summarise(
      n = dplyr::n(),
      psi_mean = mean(psi_17, na.rm = TRUE),
      psi_sd = sd(psi_17, na.rm = TRUE),
      psi_median = median(psi_17, na.rm = TRUE),
      psi_q25 = quantile(psi_17, 0.25, na.rm = TRUE),
      psi_q75 = quantile(psi_17, 0.75, na.rm = TRUE),
      kessler_mean = mean(kessler_val, na.rm = TRUE),
      wemwbs_mean = mean(wemwbs_val, na.rm = TRUE),
      age23_depression_case_rate = mean(depression_case_23, na.rm = TRUE)
    )

  covariates <- c("age17_age_z", "sex_z", "income_quintile_z", "sampling_stratum")
  covariates <- covariates[covariates %in% names(bal)]
  covariates <- covariates[vapply(bal[covariates], function(x) sum(!is.na(x)) > 100 && has_variation(x), logical(1))]
  cov_str <- if (length(covariates) > 0) paste0(" + ", paste(covariates, collapse = " + ")) else ""

  fit_fn <- function(formula_obj) {
    vars <- all.vars(formula_obj)
    dat <- bal %>%
      dplyr::filter(
        dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)),
        !is.na(GOVWT2), GOVWT2 > 0,
        !is.na(PTTYPE2), !is.na(SPTN00)
      )

    outcome_name <- vars[[1]]
    if (nrow(dat) < 100 || dplyr::n_distinct(dat[[outcome_name]]) < 2) {
      return(NULL)
    }

    design <- survey::svydesign(
      ids = ~SPTN00,
      strata = ~PTTYPE2,
      weights = ~GOVWT2,
      data = dat,
      nest = TRUE
    )

    fit <- suppressWarnings(survey::svyglm(formula_obj, design = design, family = binomial()))
    attr(fit, "balance_obs") <- dat[[outcome_name]]
    attr(fit, "balance_weights") <- dat$GOVWT2
    fit
  }

  candidate_outcomes <- c(
    depression_case_23 = paste0("Age 23 probable depressive case (PHQ-2 >= 3)", suffix)
  )

  results_list <- lapply(names(candidate_outcomes), function(outcome_var) {
    f_distress <- as.formula(paste0(outcome_var, " ~ kessler_17_z", cov_str))
    f_wellbeing <- as.formula(paste0(outcome_var, " ~ wemwbs_17_z", cov_str))
    f_balance <- as.formula(paste0(outcome_var, " ~ psi_17_z", cov_str))

    m_distress <- tryCatch(fit_fn(f_distress), error = function(e) NULL)
    m_wellbeing <- tryCatch(fit_fn(f_wellbeing), error = function(e) NULL)
    m_balance <- tryCatch(fit_fn(f_balance), error = function(e) NULL)

    dplyr::bind_rows(
      balance_extract_coef(m_distress, "kessler_17_z", "Distress only", candidate_outcomes[[outcome_var]]),
      balance_extract_coef(m_wellbeing, "wemwbs_17_z", "Wellbeing only", candidate_outcomes[[outcome_var]]),
      balance_extract_coef(m_balance, "psi_17_z", "Balance index", candidate_outcomes[[outcome_var]])
    )
  })

  comparison <- dplyr::bind_rows(results_list)
  if (nrow(comparison) == 0) {
    return(list(
      descriptives = desc,
      comparison = tibble::tibble(note = "No estimable models"),
      aic = tibble::tibble(note = "No estimable models")
    ))
  }

  comparison <- comparison %>%
    dplyr::mutate(missingness_handling = "Complete-case analysis")

  list(
    descriptives = desc,
    comparison = comparison,
    aic = balance_aic_table(comparison)
  )
}

### DT --> ================================================================
### DT --> MCS balance index
### DT --> psi_17 = log(GDCKESSL / GDWEMWBS), predicting age-23 probable-case
### DT --> depressive outcome.
### DT --> ================================================================

fit_mcs_balance_models <- function(mcs_obj) {
  df <- mcs_obj$wide
  if (is.null(df)) {
    message("MCS balance index: mcs_obj$wide is not available. Skipping.")
    return(balance_note_result("Variables not available"))
  }

  message("MCS balance index: using age-17 Kessler and WEMWBS variables from mcs_obj$wide.")

  primary_bal <- make_mcs_balance_subset(df)
  primary_results <- fit_mcs_balance_core(primary_bal)

  list(
    descriptives = primary_results$descriptives,
    comparison = primary_results$comparison,
    aic = primary_results$aic
  )
}

make_abcd_balance_subset <- function(df, threshold = 64) {
  bal <- df %>%
    dplyr::filter(
      !is.na(stress_val), !is.na(positive_affect_val),
      stress_val > 0, positive_affect_val > 0,
      !is.na(baseline_internalising_t)
    ) %>%
    dplyr::mutate(
      psi = log(stress_val / positive_affect_val),
      psi_z = zscore(psi),
      stress_z = zscore(stress_val),
      positive_affect_z = zscore(positive_affect_val),
      site = as.factor(site)
    ) %>%
    dplyr::filter(!is.na(followup_internalising_t)) %>%
    dplyr::mutate(
      probable_internalising_case = dplyr::if_else(followup_internalising_t >= threshold, 1, 0),
      baseline_noncase = dplyr::if_else(baseline_internalising_t < threshold, 1, 0)
    ) %>%
    dplyr::filter(baseline_noncase == 1)

  if (!"baseline_age_z" %in% names(bal) || all(is.na(bal$baseline_age_z))) bal$baseline_age_z <- zscore(bal$baseline_age)
  if (!"sex_z" %in% names(bal) || all(is.na(bal$sex_z))) bal$sex_z <- zscore(bal$sex)
  if (!"household_income_z" %in% names(bal) || all(is.na(bal$household_income_z))) bal$household_income_z <- zscore(bal$household_income)
  if (!"parent_education_z" %in% names(bal) || all(is.na(bal$parent_education_z))) bal$parent_education_z <- zscore(bal$parent_education)
  if (!"neighborhood_risk_z" %in% names(bal) || all(is.na(bal$neighborhood_risk_z))) bal$neighborhood_risk_z <- zscore(bal$neighborhood_risk)
  if (!"baseline_internalising_z" %in% names(bal) || all(is.na(bal$baseline_internalising_z))) bal$baseline_internalising_z <- zscore(bal$baseline_internalising_t)

  bal
}

fit_abcd_balance_core <- function(bal, outcome_var, outcome_label, denominator_model_label = "Positive affect only") {
  if (is.null(bal) || nrow(bal) == 0) {
    return(list(
      descriptives = tibble::tibble(note = "No estimable models"),
      comparison = tibble::tibble(note = "No estimable models"),
      aic = tibble::tibble(note = "No estimable models")
    ))
  }

  desc <- bal %>%
    dplyr::summarise(
      n = dplyr::n(),
      psi_mean = mean(psi, na.rm = TRUE),
      psi_sd = sd(psi, na.rm = TRUE),
      psi_median = median(psi, na.rm = TRUE),
      psi_q25 = quantile(psi, 0.25, na.rm = TRUE),
      psi_q75 = quantile(psi, 0.75, na.rm = TRUE),
      stress_mean = mean(stress_val, na.rm = TRUE),
      positive_affect_mean = mean(positive_affect_val, na.rm = TRUE),
      case_rate = mean(.data[[outcome_var]], na.rm = TRUE)
    )

  covariates <- c(
    "baseline_internalising_z", "baseline_age_z", "sex_z",
    "household_income_z", "parent_education_z", "neighborhood_risk_z"
  )
  covariates <- covariates[covariates %in% names(bal)]
  covariates <- covariates[vapply(bal[covariates], function(x) sum(!is.na(x)) > 100 && has_variation(x), logical(1))]
  if ("site" %in% names(bal) && dplyr::n_distinct(stats::na.omit(bal$site)) > 1) {
    covariates <- c(covariates, "site")
  }
  cov_str <- if (length(covariates) > 0) paste0(" + ", paste(covariates, collapse = " + ")) else ""

  fit_fn <- function(formula_obj) {
    vars <- all.vars(formula_obj)
    dat <- bal %>% dplyr::filter(dplyr::if_all(dplyr::all_of(vars), ~ !is.na(.)))
    outcome_name <- vars[[1]]
    if (nrow(dat) < 100 || dplyr::n_distinct(dat[[outcome_name]]) < 2) {
      return(NULL)
    }
    fit <- glm(formula_obj, data = dat, family = binomial())
    attr(fit, "balance_obs") <- dat[[outcome_name]]
    attr(fit, "balance_weights") <- rep(1, nrow(dat))
    fit
  }

  m_stress <- tryCatch(fit_fn(as.formula(paste0(outcome_var, " ~ stress_z", cov_str))), error = function(e) NULL)
  m_positive_affect <- tryCatch(fit_fn(as.formula(paste0(outcome_var, " ~ positive_affect_z", cov_str))), error = function(e) NULL)
  m_balance <- tryCatch(fit_fn(as.formula(paste0(outcome_var, " ~ psi_z", cov_str))), error = function(e) NULL)

  comparison <- dplyr::bind_rows(
    balance_extract_coef(m_stress, "stress_z", "Stress only", outcome_label),
    balance_extract_coef(m_positive_affect, "positive_affect_z", denominator_model_label, outcome_label),
    balance_extract_coef(m_balance, "psi_z", "Balance index", outcome_label)
  )

  if (nrow(comparison) == 0) {
    return(list(
      descriptives = desc,
      comparison = tibble::tibble(note = "No estimable models"),
      aic = tibble::tibble(note = "No estimable models")
    ))
  }

  comparison <- comparison %>%
    dplyr::mutate(missingness_handling = "Complete-case analysis")

  list(
    descriptives = desc,
    comparison = comparison,
    aic = balance_aic_table(comparison)
  )
}

make_abcd_balance_base <- function(
  abcd_obj,
  baseline_session,
  followup_session,
  stress_var = "cbcl_stress_sum",
  denominator_var = "positive_affect_sum"
) {
  long <- abcd_obj$long
  if (is.null(long) || nrow(long) == 0) {
    return(NULL)
  }

  required_vars <- c(
    "participant_id", "session_id", "age", "sex", "site",
    "household_income", "parent_education", "neighborhood_risk",
    "internalising_t", stress_var, denominator_var
  )
  if (!all(required_vars %in% names(long))) {
    return(NULL)
  }

  baseline_df <- long %>%
    dplyr::filter(session_id == baseline_session) %>%
    dplyr::transmute(
      participant_id,
      baseline_age = age,
      sex = sex,
      site = site,
      household_income = household_income,
      parent_education = parent_education,
      neighborhood_risk = neighborhood_risk,
      baseline_internalising_t = internalising_t,
      stress_val = .data[[stress_var]],
      positive_affect_val = .data[[denominator_var]]
    ) %>%
    dplyr::distinct(participant_id, .keep_all = TRUE)

  followup_df <- long %>%
    dplyr::filter(session_id == followup_session) %>%
    dplyr::transmute(
      participant_id,
      followup_internalising_t = internalising_t
    ) %>%
    dplyr::distinct(participant_id, .keep_all = TRUE)

  baseline_df %>%
    dplyr::left_join(followup_df, by = "participant_id") %>%
    dplyr::mutate(
      baseline_age_z = zscore(baseline_age),
      sex_z = zscore(sex),
      household_income_z = zscore(household_income),
      parent_education_z = zscore(parent_education),
      neighborhood_risk_z = zscore(neighborhood_risk),
      baseline_internalising_z = zscore(baseline_internalising_t),
      site = as.factor(site)
    )
}

### DT --> ================================================================
### DT --> ABCD balance index
### DT --> psi = log(CBCL stress / Positive Affect), using ses-03A stress
### DT --> and ses-03A youth-reported positive affect to predict later
### DT --> internalising caseness at ses-05A.
### DT --> ================================================================

fit_abcd_balance_models <- function(abcd_obj, config) {
  outcome_label <- "Age 15 probable internalising case (CBCL >= 64; age-13 baseline < 64)"
  bal_base <- make_abcd_balance_base(
    abcd_obj = abcd_obj,
    baseline_session = "ses-03A",
    followup_session = config$analysis$abcd_followup_session,
    stress_var = "cbcl_stress_sum",
    denominator_var = "positive_affect_sum"
  )

  if (is.null(bal_base) || nrow(bal_base) == 0) {
    message("ABCD balance index: required age-13 variables not found. Skipping.")
    return(balance_note_result("Variables not available"))
  }

  results <- fit_abcd_balance_core(
    make_abcd_balance_subset(bal_base, threshold = 64),
    outcome_var = "probable_internalising_case",
    outcome_label = outcome_label,
    denominator_model_label = "Positive affect only"
  )

  if (!"note" %in% names(results$descriptives)) {
    results$descriptives <- results$descriptives %>%
      dplyr::mutate(outcome = outcome_label, .before = 1)
  }

  results
}
