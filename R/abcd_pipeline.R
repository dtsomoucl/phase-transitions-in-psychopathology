### DT --> ABCD-specific mappings, construct derivation, longitudinal summaries, and lagged predictive dataset assembly.
abcd_ethnicity_5cat <- function(hispanic, white, black, asian_count, total_race_count) {
  factor(
    dplyr::case_when(
      hispanic == 1 ~ "Hispanic/Latino",
      total_race_count == 1 & white == 1 ~ "White",
      total_race_count == 1 & black == 1 ~ "Black",
      asian_count >= 1 & asian_count == total_race_count ~ "Asian",
      total_race_count >= 1 ~ "Other/Multiracial",
      TRUE ~ NA_character_
    ),
    levels = c("White", "Black", "Hispanic/Latino", "Asian", "Other/Multiracial")
  )
}

abcd_variable_mapping <- function(config) {
  tibble(
    cohort = "ABCD",
    category = c(
      "id_wave", "id_wave", rep("covariate", 6),
      rep("outcome", 6),
      "field_like", "field_like", "field_like",
      "field_sensitivity", "field_sensitivity", "field_sensitivity",
      "precision_like", "precision_sensitivity", "precision_sensitivity", "precision_sensitivity",
      "context", "context", "context", "context", "context",
      "prs"
    ),
    construct = c(
      "participant_id", "session_id", "age", "sex", "site", "household_income", "parent_education", "neighborhood_disadvantage",
      "internalising", "withdrawn_depressed", "depressive_problems", "youth_bpm_internalising", "youth_bpm_attention", "youth_depression_impairment",
      "school_environment", "school_involvement", "school_disengagement",
      "bas_drive", "bas_fun_seeking", "bas_reward_responsiveness",
      "flanker_control", "pattern_comparison_processing_speed", "crystallized_composite",
      "card_sort_flexibility_legacy", "list_sorting_working_memory_legacy",
      "family_cohesion", "family_conflict", "sleep_duration", "social_jetlag", "neighborhood_adi_national_percentile",
      "prs_local_status"
    ),
    file = fs::path(
      config$raw_data$abcd_dir,
      c(
        "mh_p_cbcl.tsv", "mh_p_cbcl.tsv", "mh_p_cbcl.tsv", "ab_p_demo.tsv", "ab_p_demo.tsv",
        "ab_g_stc.tsv", "ab_g_dyn.tsv", "le_l_adi.tsv",
        "mh_p_cbcl.tsv", "mh_p_cbcl.tsv", "mh_p_cbcl.tsv", "mh_y_bpm.tsv", "mh_y_bpm.tsv", "mh_y_ksads__dep.parquet",
        "fc_y_srpf.tsv", "fc_y_srpf.tsv", "fc_y_srpf.tsv",
        "mh_y_bisbas.tsv", "mh_y_bisbas.tsv", "mh_y_bisbas.tsv",
        "nc_y_nihtb.tsv", "nc_y_nihtb.tsv", "nc_y_nihtb.tsv", "nc_y_nihtb.tsv", "nc_y_nihtb.tsv",
        "fc_y_fes.tsv", "fc_y_fes.tsv", "ph_y_mctq.tsv", "ph_y_mctq.tsv", "le_l_adi.tsv",
        NA_character_
      )
    ),
    variable = c(
      "participant_id", "session_id", "mh_p_cbcl_age", "ab_g_stc__cohort_sex", "ab_g_dyn__design_site", "ab_p_demo__income__hhold_001", "ab_p_demo__edu__slf_001__v02", "le_l_adi__addr1__national_prcnt",
      "mh_p_cbcl__synd__int_tscore", "mh_p_cbcl__synd__wthdep_tscore", "mh_p_cbcl__dsm__dep_tscore", "mh_y_bpm__int_tscore", "mh_y_bpm__attn_tscore", "mh_y_ksads__dep__impfunct__pres_sx",
      "fc_y_srpf__env_mean", "fc_y_srpf__involv_mean", "fc_y_srpf__dis_mean",
      "mh_y_bisbas__bas__dr_sum", "mh_y_bisbas__bas__fs_sum", "mh_y_bisbas__bas__rr_sum",
      "nc_y_nihtb__flnkr__fullcor_tscore",
      "nc_y_nihtb__pttcp__fullcor_tscore", "nc_y_nihtb__comp__cryst__fullcorr_tscore",
      "nc_y_nihtb__crdst__fullcorr_tscore", "nc_y_nihtb__lswmt__fullcor_tscore",
      "fc_y_fes__cohes_mean", "fc_y_fes__confl_mean", "ph_y_mctq__sleep_dur", "ph_y_mctq__socjl_absl", "le_l_adi__addr1__national_prcnt",
      "No local ABCD PRS/PGI file located in this workspace"
    ),
    notes = c(
      "ABCD participant identifier", "ABCD session identifier", "Table-specific age", "Stable cohort sex", "Assessment site", "Parent-reported household income band", "Parent education version 2 when available", "Linked area deprivation index percentile (address 1)",
      "Primary outcome for ABCD longitudinal analyses", "Withdrawal/disengagement proxy", "Secondary depressive outcome", "Youth-reported internalising sensitivity outcome", "Youth-reported attention-problems sensitivity outcome", "Youth KSADS depression-related functional impairment; used in the secondary balance-index analysis at ses-04A",
      "School climate/protection", "School involvement", "Reverse-scored in the primary field composite",
      "Alternative field-like candidate", "Alternative field-like candidate", "Alternative field-like candidate",
      "Primary executive-control proxy used in the main models",
      "Precision-proxy sensitivity: Pattern Comparison Processing Speed (good ses-02A coverage)",
      "Precision-proxy sensitivity: Crystallized composite (good ses-02A coverage)",
      "Replaced by pttcp/crystallized at ses-02A: only 5 valid cases", "Replaced by pttcp/crystallized at ses-02A: only 14 valid cases",
      "Protective family context", "Risk family context", "Higher values are protective", "Higher values are risk", "Higher values indicate greater area deprivation",
      "Local PRS assets were not detected for ABCD"
    )
  )
}

build_abcd_dataset <- function(config) {
  abcd_dir <- config$raw_data$abcd_dir

  demo <- read_abcd_tsv_any(
    fs::path(abcd_dir, "ab_p_demo.tsv"),
    c(
      "participant_id", "session_id", "ab_p_demo__income__hhold_001",
      "ab_p_demo__income__hhold_001__v01",
      "ab_p_demo__edu__slf_001", "ab_p_demo__edu__slf_001__v01", "ab_p_demo__edu__slf_001__v02",
      "ab_p_demo__ethn_001", "ab_p_demo__ethn__slf_001",
      "ab_p_demo__race_001___10", "ab_p_demo__race_001___11", "ab_p_demo__race_001___12",
      "ab_p_demo__race_001___13", "ab_p_demo__race_001___14", "ab_p_demo__race_001___15",
      "ab_p_demo__race_001___16", "ab_p_demo__race_001___17", "ab_p_demo__race_001___18",
      "ab_p_demo__race_001___19", "ab_p_demo__race_001___20", "ab_p_demo__race_001___21",
      "ab_p_demo__race_001___22", "ab_p_demo__race_001___23", "ab_p_demo__race_001___24",
      "ab_p_demo__race_001___25"
    )
  )

  cohort_static <- read_abcd_tsv_any(
    fs::path(abcd_dir, "ab_g_stc.tsv"),
    c("participant_id", "ab_g_stc__cohort_sex")
  )

  visit_design <- read_abcd_tsv_any(
    fs::path(abcd_dir, "ab_g_dyn.tsv"),
    c("participant_id", "session_id", "ab_g_dyn__design_site")
  )

  cbcl <- read_abcd_tsv_any(
    fs::path(abcd_dir, "mh_p_cbcl.tsv"),
    c(
      "participant_id", "session_id", "mh_p_cbcl_age",
      "mh_p_cbcl__strs_sum",
      "mh_p_cbcl__synd__int_tscore", "mh_p_cbcl__synd__wthdep_tscore",
      "mh_p_cbcl__dsm__dep_tscore", "mh_p_cbcl__synd__attn_tscore"
    )
  )

  bpm_y <- read_abcd_tsv_any(
    fs::path(abcd_dir, "mh_y_bpm.tsv"),
    c(
      "participant_id", "session_id", "mh_y_bpm_age",
      "mh_y_bpm__int_tscore", "mh_y_bpm__attn_tscore"
    )
  )

  ksads_dep_subset_path <- fs::path(abcd_dir, "mh_y_ksads__dep_impair_subset.tsv")
  ksads_dep <- if (file.exists(ksads_dep_subset_path)) {
    read_abcd_tsv_any(
      ksads_dep_subset_path,
      c(
        "participant_id", "session_id",
        "mh_y_ksads__dep__impfunct__pres_sx"
      )
    )
  } else {
    read_abcd_parquet_any(
      fs::path(abcd_dir, "mh_y_ksads__dep.parquet"),
      c(
        "participant_id", "session_id",
        "mh_y_ksads__dep__impfunct__pres_sx"
      )
    )
  }

  bisbas <- read_abcd_tsv_any(
    fs::path(abcd_dir, "mh_y_bisbas.tsv"),
    c(
      "participant_id", "session_id", "mh_y_bisbas_age",
      "mh_y_bisbas__bas__dr_sum", "mh_y_bisbas__bas__fs_sum",
      "mh_y_bisbas__bas__rr_sum", "mh_y_bisbas__bas__rr_sum__v01",
      "mh_y_bisbas__bis_sum", "mh_y_bisbas__bis_sum__v01"
    )
  )

  school <- read_abcd_tsv_any(
    fs::path(abcd_dir, "fc_y_srpf.tsv"),
    c(
      "participant_id", "session_id", "fc_y_srpf_age",
      "fc_y_srpf__dis_mean", "fc_y_srpf__env_mean", "fc_y_srpf__involv_mean"
    )
  )

  nihtb <- read_abcd_tsv_any(
    fs::path(abcd_dir, "nc_y_nihtb.tsv"),
    c(
      "participant_id", "session_id",
      "nc_y_nihtb__comp__fluid__fullcorr_tscore",
      "nc_y_nihtb__crdst__fullcorr_tscore",
      "nc_y_nihtb__flnkr__fullcor_tscore",
      "nc_y_nihtb__lswmt__fullcor_tscore",
      ### DT ---> Precision-proxy sensitivity replacements (good ses-02A coverage)
      "nc_y_nihtb__pttcp__fullcor_tscore",           ### DT ---> Pattern Comparison Processing Speed
      "nc_y_nihtb__comp__cryst__fullcorr_tscore"     ### DT ---> Crystallized composite
    )
  )

  bdefs <- read_abcd_tsv_any(
    fs::path(abcd_dir, "nc_p_bdefs.tsv"),
    c("participant_id", "session_id", "nc_p_bdefs_sum")
  )

  fes <- read_abcd_tsv_any(
    fs::path(abcd_dir, "fc_y_fes.tsv"),
    c(
      "participant_id", "session_id", "fc_y_fes__cohes_mean", "fc_y_fes__confl_mean"
    )
  )

  sleep <- read_abcd_tsv_any(
    fs::path(abcd_dir, "ph_y_mctq.tsv"),
    c(
      "participant_id", "session_id", "ph_y_mctq__sleep_dur", "ph_y_mctq__socjl_absl"
    )
  )

  cyber <- read_abcd_tsv_any(
    fs::path(abcd_dir, "mh_y_cb.tsv"),
    c(
      "participant_id", "session_id", "mh_y_cb_001a", "mh_y_cb_001a__01", "mh_y_cb_001a__01__01"
    )
  )

  adi <- read_abcd_tsv_any(
    fs::path(abcd_dir, "le_l_adi.tsv"),
    c(
      "participant_id", "session_id", "le_l_adi__addr1__national_prcnt"
    )
  ) %>%
    select(-session_id) %>%
    collapse_by_id("participant_id")

  long <- cbcl %>%
    left_join(demo, by = c("participant_id", "session_id")) %>%
    left_join(cohort_static, by = "participant_id") %>%
    left_join(visit_design, by = c("participant_id", "session_id")) %>%
    left_join(bpm_y, by = c("participant_id", "session_id")) %>%
    left_join(ksads_dep, by = c("participant_id", "session_id")) %>%
    left_join(bisbas, by = c("participant_id", "session_id")) %>%
    left_join(school, by = c("participant_id", "session_id")) %>%
    left_join(nihtb, by = c("participant_id", "session_id")) %>%
    left_join(bdefs, by = c("participant_id", "session_id")) %>%
    left_join(fes, by = c("participant_id", "session_id")) %>%
    left_join(sleep, by = c("participant_id", "session_id")) %>%
    left_join(cyber, by = c("participant_id", "session_id")) %>%
    left_join(adi, by = "participant_id") %>%
    mutate(
      wave_year = session_to_wave_year(session_id),
      age = clean_numeric(mh_p_cbcl_age),
      sex = clean_numeric(ab_g_stc__cohort_sex),
      site = as.factor(clean_numeric(ab_g_dyn__design_site)),
      household_income = first_non_missing(
        clean_numeric(ab_p_demo__income__hhold_001),
        clean_numeric(ab_p_demo__income__hhold_001__v01)
      ),
      parent_education = first_non_missing(
        clean_numeric(ab_p_demo__edu__slf_001__v02),
        clean_numeric(ab_p_demo__edu__slf_001__v01),
        clean_numeric(ab_p_demo__edu__slf_001)
      ),
      hispanic_ethnicity = first_non_missing(
        clean_numeric(ab_p_demo__ethn_001),
        clean_numeric(ab_p_demo__ethn__slf_001)
      ),
      race_white = clean_numeric(ab_p_demo__race_001___10),
      race_black = clean_numeric(ab_p_demo__race_001___11),
      race_aian = clean_numeric(ab_p_demo__race_001___12),
      race_asian_indian = clean_numeric(ab_p_demo__race_001___18),
      race_asian_chinese = clean_numeric(ab_p_demo__race_001___19),
      race_asian_filipino = clean_numeric(ab_p_demo__race_001___20),
      race_asian_japanese = clean_numeric(ab_p_demo__race_001___21),
      race_asian_korean = clean_numeric(ab_p_demo__race_001___22),
      race_asian_vietnamese = clean_numeric(ab_p_demo__race_001___23),
      race_asian_other = clean_numeric(ab_p_demo__race_001___24),
      race_other = clean_numeric(ab_p_demo__race_001___25),
      race_total_count = rowSums(across(c(
        race_white, race_black, race_aian,
        ab_p_demo__race_001___13, ab_p_demo__race_001___14, ab_p_demo__race_001___15,
        ab_p_demo__race_001___16, ab_p_demo__race_001___17,
        race_asian_indian, race_asian_chinese, race_asian_filipino, race_asian_japanese,
        race_asian_korean, race_asian_vietnamese, race_asian_other, race_other
      ), clean_numeric), na.rm = TRUE),
      race_asian_count = rowSums(across(c(
        race_asian_indian, race_asian_chinese, race_asian_filipino, race_asian_japanese,
        race_asian_korean, race_asian_vietnamese, race_asian_other
      )), na.rm = TRUE),
      ethnicity_5cat = abcd_ethnicity_5cat(
        hispanic = hispanic_ethnicity,
        white = race_white,
        black = race_black,
        asian_count = race_asian_count,
        total_race_count = race_total_count
      ),
      internalising_t = clean_numeric(mh_p_cbcl__synd__int_tscore),
      cbcl_stress_sum = clean_numeric(mh_p_cbcl__strs_sum),
      withdrawn_t = clean_numeric(mh_p_cbcl__synd__wthdep_tscore),
      depressive_t = clean_numeric(mh_p_cbcl__dsm__dep_tscore),
      bpm_internalising_t = clean_numeric(mh_y_bpm__int_tscore),
      bpm_attention_t = clean_numeric(mh_y_bpm__attn_tscore),
      depression_impair_raw = suppressWarnings(as.numeric(mh_y_ksads__dep__impfunct__pres_sx)),
      depression_impair_case = case_when(
        depression_impair_raw == 1 ~ 1,
        depression_impair_raw %in% c(0, 888) ~ 0,
        TRUE ~ NA_real_
      ),
      bas_drive = clean_numeric(mh_y_bisbas__bas__dr_sum),
      bas_fun = clean_numeric(mh_y_bisbas__bas__fs_sum),
      bas_reward = first_non_missing(
        clean_numeric(mh_y_bisbas__bas__rr_sum),
        clean_numeric(mh_y_bisbas__bas__rr_sum__v01)
      ),
      school_env = clean_numeric(fc_y_srpf__env_mean),
      school_involvement = clean_numeric(fc_y_srpf__involv_mean),
      school_disengagement = clean_numeric(fc_y_srpf__dis_mean),
      nihtb_flanker = clean_numeric(nc_y_nihtb__flnkr__fullcor_tscore),
      nihtb_cardsort = clean_numeric(nc_y_nihtb__crdst__fullcorr_tscore),
      nihtb_workmem = clean_numeric(nc_y_nihtb__lswmt__fullcor_tscore),
      ### DT ---> Precision-proxy sensitivity replacements
      nihtb_processing_speed = clean_numeric(nc_y_nihtb__pttcp__fullcor_tscore),
      nihtb_crystallized = clean_numeric(nc_y_nihtb__comp__cryst__fullcorr_tscore),
      nihtb_fluid = clean_numeric(nc_y_nihtb__comp__fluid__fullcorr_tscore),
      bdefs_total = clean_numeric(nc_p_bdefs_sum),
      family_cohesion = clean_numeric(fc_y_fes__cohes_mean),
      family_conflict = clean_numeric(fc_y_fes__confl_mean),
      sleep_duration = clean_numeric(ph_y_mctq__sleep_dur),
      social_jetlag = clean_numeric(ph_y_mctq__socjl_absl),
      cyberbullied_12m = first_non_missing(
        clean_numeric(mh_y_cb_001a__01),
        if_else(clean_numeric(mh_y_cb_001a__01__01) > 0, 1, 0, missing = NA_real_)
      ),
      neighborhood_adi = clean_numeric(le_l_adi__addr1__national_prcnt),
      reverse_school_disengagement = -clean_numeric(fc_y_srpf__dis_mean),
      reverse_bdefs_total = -clean_numeric(nc_p_bdefs_sum)
    ) %>%
    group_by(session_id) %>%
    mutate(
      field_index = row_mean_standardised(
        pick(everything()),
        c("school_env", "school_involvement", "reverse_school_disengagement")
      ),
      precision_index = zscore(nihtb_flanker),
      precision_cardsort_index = zscore(nihtb_cardsort),
      precision_workmem_index = zscore(nihtb_workmem),
      ### DT ---> Precision-proxy sensitivity replacements
      precision_processing_speed_index = zscore(nihtb_processing_speed),
      precision_crystallized_index = zscore(nihtb_crystallized),
      family_context = zscore(family_cohesion) - zscore(family_conflict),
      sleep_risk = zscore(-sleep_duration) + zscore(social_jetlag),
      peer_adversity = zscore(cyberbullied_12m),
      neighborhood_risk = zscore(neighborhood_adi),
      internalising_z = zscore(internalising_t),
      withdrawn_z = zscore(withdrawn_t),
      bpm_internalising_z = zscore(bpm_internalising_t),
      bpm_attention_z = zscore(bpm_attention_t)
    ) %>%
    ungroup() %>%
    group_by(participant_id) %>%
    tidyr::fill(ethnicity_5cat, .direction = "downup") %>%
    ungroup()

  session_pairs <- expand.grid(
    baseline_session = sort(unique(long$session_id)),
    followup_session = sort(unique(long$session_id)),
    stringsAsFactors = FALSE
  ) %>%
    filter(baseline_session < followup_session) %>%
    mutate(
      complete_cases = map2_int(baseline_session, followup_session, function(b, f) {
        pair_df <- long %>%
          filter(session_id %in% c(b, f)) %>%
          select(participant_id, session_id, internalising_t, field_index, precision_index) %>%
          pivot_wider(
            id_cols = participant_id,
            names_from = session_id,
            values_from = c(internalising_t, field_index, precision_index),
            names_sep = "__"
          )
        sum(
          complete.cases(
            pair_df[[paste0("internalising_t__", b)]],
            pair_df[[paste0("internalising_t__", f)]],
            pair_df[[paste0("field_index__", b)]],
            pair_df[[paste0("precision_index__", b)]]
          )
        )
      })
    ) %>%
    arrange(desc(complete_cases), baseline_session, followup_session)

  requested_pair <- session_pairs %>%
    filter(
      baseline_session == config$analysis$abcd_baseline_session,
      followup_session == config$analysis$abcd_followup_session
    )

  chosen_pair <- if (nrow(requested_pair) == 1 && requested_pair$complete_cases > 0) {
    requested_pair
  } else {
    session_pairs %>% filter(complete_cases > 0) %>% slice(1)
  }

  baseline_session <- chosen_pair$baseline_session[[1]]
  followup_session <- chosen_pair$followup_session[[1]]

  predictive <- long %>%
    filter(session_id %in% c(baseline_session, followup_session)) %>%
    select(
      participant_id, session_id, age, ethnicity_5cat, household_income, parent_education,
      sex, site, neighborhood_adi, neighborhood_risk,
      internalising_t, internalising_z, cbcl_stress_sum, withdrawn_t, withdrawn_z,
      bas_drive,
      bpm_internalising_t, bpm_internalising_z, bpm_attention_t, bpm_attention_z,
      field_index, precision_index,
      precision_cardsort_index, precision_workmem_index,
      precision_processing_speed_index, precision_crystallized_index,
      family_context, sleep_risk,
      peer_adversity
    ) %>%
    pivot_wider(
      id_cols = participant_id,
      names_from = session_id,
      values_from = -participant_id,
      names_sep = "__"
    ) %>%
    transmute(
      participant_id = participant_id,
      baseline_age = .data[[paste0("age__", baseline_session)]],
      ethnicity_5cat = .data[[paste0("ethnicity_5cat__", baseline_session)]],
      sex = .data[[paste0("sex__", baseline_session)]],
      site = .data[[paste0("site__", baseline_session)]],
      household_income = .data[[paste0("household_income__", baseline_session)]],
      parent_education = .data[[paste0("parent_education__", baseline_session)]],
      neighborhood_adi = .data[[paste0("neighborhood_adi__", baseline_session)]],
      neighborhood_risk = .data[[paste0("neighborhood_risk__", baseline_session)]],
      baseline_internalising_t = .data[[paste0("internalising_t__", baseline_session)]],
      baseline_internalising_z = .data[[paste0("internalising_z__", baseline_session)]],
      baseline_stress_sum = .data[[paste0("cbcl_stress_sum__", baseline_session)]],
      baseline_bas_drive = .data[[paste0("bas_drive__", baseline_session)]],
      followup_internalising_t = .data[[paste0("internalising_t__", followup_session)]],
      followup_internalising_z = .data[[paste0("internalising_z__", followup_session)]],
      baseline_bpm_internalising_z = .data[[paste0("bpm_internalising_z__", baseline_session)]],
      followup_bpm_internalising_z = .data[[paste0("bpm_internalising_z__", followup_session)]],
      baseline_bpm_attention_z = .data[[paste0("bpm_attention_z__", baseline_session)]],
      followup_bpm_attention_z = .data[[paste0("bpm_attention_z__", followup_session)]],
      followup_withdrawn_t = .data[[paste0("withdrawn_t__", followup_session)]],
      followup_withdrawn_z = .data[[paste0("withdrawn_z__", followup_session)]],
      field_index = .data[[paste0("field_index__", baseline_session)]],
      precision_index = .data[[paste0("precision_index__", baseline_session)]],
      precision_cardsort_index = .data[[paste0("precision_cardsort_index__", baseline_session)]],
      precision_workmem_index = .data[[paste0("precision_workmem_index__", baseline_session)]],
      ### DT ---> Precision-proxy sensitivity replacements
      precision_processing_speed_index = .data[[paste0("precision_processing_speed_index__", baseline_session)]],
      precision_crystallized_index = .data[[paste0("precision_crystallized_index__", baseline_session)]],
      family_context = .data[[paste0("family_context__", baseline_session)]],
      sleep_risk = .data[[paste0("sleep_risk__", baseline_session)]],
      peer_adversity = .data[[paste0("peer_adversity__", baseline_session)]]
    )

  wide <- predictive %>%
    mutate(
      cbcl_stress_baseline = baseline_stress_sum,
      bas_drive_baseline = baseline_bas_drive,
      internalising_followup = followup_internalising_z,
      baseline_age_z = zscore(baseline_age),
      sex_z = zscore(sex),
      household_income_z = zscore(household_income),
      parent_education_z = zscore(parent_education),
      neighborhood_risk_z = zscore(neighborhood_risk),
      site = as.factor(site)
    )

  dep_impair_ses04 <- long %>%
    filter(session_id == "ses-04A") %>%
    select(participant_id, depression_impair_case) %>%
    distinct(participant_id, .keep_all = TRUE) %>%
    rename(depression_impair_case_ses04 = depression_impair_case)

  wide <- wide %>%
    left_join(dep_impair_ses04, by = "participant_id")

  reliability_inputs <- list(
    field = long %>%
      filter(session_id == baseline_session) %>%
      select(
        participant_id,
        school_env,
        school_involvement,
        reverse_school_disengagement
      ),
    precision = long %>%
      filter(session_id == baseline_session) %>%
      select(
        participant_id,
        nihtb_flanker
      )
  )

  list(
    long = long,
    predictive = predictive,
    wide = wide,
    baseline_session = baseline_session,
    followup_session = followup_session,
    session_pairs = session_pairs,
    reliability_inputs = reliability_inputs
  )
}
