### DT --> Reliability summaries for the theory-motivated composite scales.

compute_standardised_reliability <- function(data, id_var, items, cohort, wave, composite, exclude_values = NULL) {
  item_df <- data %>%
    select(all_of(c(id_var, items))) %>%
    distinct()

  complete_items <- item_df %>%
    select(all_of(items)) %>%
    mutate(across(everything(), clean_numeric)) %>%
    mutate(across(everything(), ~ if (!is.null(exclude_values)) replace(., . %in% exclude_values, NA_real_) else .)) %>%
    tidyr::drop_na()

  out <- tibble(
    cohort = cohort,
    wave = wave,
    composite = composite,
    n_total = nrow(item_df),
    n_complete = nrow(complete_items),
    item_count = length(items),
    items = paste(items, collapse = "; "),
    alpha_std = NA_real_,
    omega_total = NA_real_
  )

  if (nrow(complete_items) < 2 || length(items) < 2) {
    return(out)
  }

  alpha_fit <- suppressWarnings(psych::alpha(complete_items, check.keys = FALSE, warnings = FALSE))
  omega_fit <- suppressMessages(suppressWarnings(
    psych::omega(complete_items, nfactors = 1, plot = FALSE, warnings = FALSE)
  ))

  out %>%
    mutate(
      alpha_std = unname(alpha_fit$total$std.alpha %||% NA_real_),
      omega_total = unname(omega_fit$omega.tot %||% NA_real_)
    )
}

build_mcs_composite_reliability <- function(mcs_obj) {
  inputs <- mcs_obj$reliability_inputs
  bind_rows(
    compute_standardised_reliability(
      data = inputs$field,
      id_var = "MCSID",
      items = c(
        "school_happy_pos",
        "schoolwork_happy_pos",
        "friends_happy_pos",
        "try_best_pos",
        "school_interest_pos",
        "school_waste_reverse_pos"
      ),
      cohort = "MCS",
      wave = "age14",
      composite = "Motivation/engagement composite"
    ),
    compute_standardised_reliability(
      data = inputs$precision,
      id_var = "MCSID",
      items = c("quality_decision_pos"),
      cohort = "MCS",
      wave = "age14",
      composite = "Executive-control primary proxy"
    )
  )
}

build_abcd_composite_reliability <- function(abcd_obj) {
  inputs <- abcd_obj$reliability_inputs
  bind_rows(
    compute_standardised_reliability(
      data = inputs$field,
      id_var = "participant_id",
      items = c(
        "school_env",
        "school_involvement",
        "reverse_school_disengagement"
      ),
      cohort = "ABCD",
      wave = abcd_obj$baseline_session,
      composite = "Motivation/engagement composite"
    ),
    compute_standardised_reliability(
      data = inputs$precision,
      id_var = "participant_id",
      items = c("nihtb_flanker"),
      cohort = "ABCD",
      wave = abcd_obj$baseline_session,
      composite = "Cognitive-control primary proxy"
    ),
    compute_standardised_reliability(
      data = inputs$positive_affect,
      id_var = "participant_id",
      items = sprintf("mh_y_pai_%03d", 1:9),
      cohort = "ABCD",
      wave = "ses-03A",
      composite = "Positive-affect balance denominator",
      exclude_values = c(777, 999)
    )
  )
}
