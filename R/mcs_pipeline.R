### DT --> MCS-specific mappings, harmonised distress proxy construction, age-14 baseline predictors, and PRS assembly.
mcs_pttype2_factor <- function(x) {
  factor(
    x,
    levels = c(1, 2, 3, 4, 5, 6, 7, 8, 9),
    labels = c(
      "England - Advantaged",
      "England - Disadvantaged",
      "England - Ethnic",
      "Wales - Advantaged",
      "Wales - Disadvantaged",
      "Scotland - Advantaged",
      "Scotland - Disadvantaged",
      "Northern Ireland - Advantaged",
      "Northern Ireland - Disadvantaged"
    )
  )
}

mcs_ethnicity_6cat <- function(x) {
  factor(
    dplyr::case_when(
      x == 1 ~ "White",
      x == 2 ~ "Mixed",
      x == 3 ~ "Indian",
      x == 4 ~ "Pakistani/Bangladeshi",
      x == 5 ~ "Black/Black British",
      x == 6 ~ "Other",
      TRUE ~ NA_character_
    ),
    levels = c("White", "Mixed", "Indian", "Pakistani/Bangladeshi", "Black/Black British", "Other")
  )
}

mcs_variable_mapping <- function(config) {
  pc_constructs <- paste0("ancestry_pc", 1:20)
  pc_variables <- paste0("pc", 1:20)

  tibble(
    cohort = "MCS",
    category = c(
      "id_wave", "id_wave",
      rep("covariate", 28),
      rep("outcome", 10),
      rep("field_like", 6),
      rep("precision_like", 1),
      rep("precision_sensitivity", 2),
      rep("context", 4),
      rep("prs", 12)
    ),
    construct = c(
      "MCSID", "wave",
      "sex_age14", "age_age14", "pttype2_age14_stratum", "oecd_income_quintile_age14", "mcs_design_stratum", "mcs_design_psu", "mcs6_weight_uk", "mcs7_weight_uk", pc_constructs,
      "age11_internalising_proxy", "age14_internalising_proxy", "age17_psychological_distress", "age17_wellbeing",
      "age14_mfq_depressive_symptoms", "age23_psychological_distress", "age23_anxiety", "age23_depression", "age23_loneliness", "age23_wellbeing",
      "school_happiness_age14", "schoolwork_happiness_age14", "friends_happiness_age14", "try_best_at_school_age14", "school_interest_age14", "reverse_school_waste_age14",
      "quality_of_decision_age14",
      "delay_aversion_age14", "word_activity_score_age14",
      "family_conflict_age14", "sleep_latency_age14", "peer_school_trouble_age14", "social_support_age14",
      "depression_prs", "mdd_prs", "anxiety_prs", "adhd_prs", "externalising_prs", "cognition_pgi", "asd_pgi",
      "agreeableness_pgi", "conscientiousness_pgi", "extraversion_pgi", "neuroticism_pgi", "openness_pgi"
    ),
    file = c(
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      NA_character_,
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_family_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      rep(fs::path(config$raw_data$mcs_dir, "cls_mcs_pgi_final_v1_1_mcsid.sav"), 20),
      fs::path(config$raw_data$mcs_dir, "mcs5_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs7_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs7_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      rep(fs::path(config$raw_data$mcs_dir, "mcs8_23y_cm_derived.sav"), 5),
      rep(fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"), 6),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_cognitive_assessment.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      rep(fs::path(config$raw_data$mcs_dir, "cls_mcs_pgi_final_v1_1_mcsid.sav"), 12)
    ),
    variable = c(
      "MCSID", "age11/age14/age17/age23",
      "FCCSEX00", "FCCAGE00", "PTTYPE2", "FOECDUK0", "PTTYPE2", "SPTN00", "FOVWT2", "GOVWT2", pc_variables,
      "EEMOTION + EPEER", "FEMOTION + FPEER", "GDCKESSL", "GDWEMWBS", "sum(FCMDSA00:FCMDSM00)", "HDKESSLER6", "HDGAD2", "HDPHQ2", "HDLONELINESS", "HDWEMWB",
      "FCSCHL00", "FCSCWK00", "FCFRNS00", "FCSCBE00", "FCSINT00", "FCSCWA00",
      "FCGTQOFDM",
      "FCGTDELAY", "FCWRDSC",
      "FCQUAM00 + FCQUAF00", "FCSLLN00", "FCPETR00", "FCSAFF00",
      "mcs_depression_p1_v1", "mcs_mdd_p1_v1", "mcs_anxiety_p1_v1", "mcs_adhd_p1_v1", "mcs_externalising_p5e08_v1", "mcs_cognition_p1_v1", "mcs_asd_p1_v1",
      "mcs_agreeableness_p1_v1", "mcs_consc_p1_v1", "mcs_extraversion_p1_v1", "mcs_neuroticism_p1_v1", "mcs_openness_p1_v1"
    ),
    notes = c(
      "Family identifier used throughout analytic joins", "Wave labels assigned in pipeline",
      "Sex at age 14 baseline", "Age at age 14 sweep", "Age-14 PTTYPE2 sampling stratum used as a sociodemographic covariate", "OECD equivalised income quintile at age 14",
      "Survey stratum for design-based inference", "Fieldwork point / PSU for design-based inference", "MCS6 whole-UK weight used when age 14 is the exposure wave", "MCS7 whole-UK weight used when age 17 is the exposure wave",
      rep("Genetic ancestry principal component used as a covariate in genetic-score models", 20),
      "Parent-reported internalising proxy at age 11", "Parent-reported internalising proxy at age 14", "Age 17 standalone psychological distress outcome", "Age 17 standalone wellbeing outcome retained separately", "Age 14 short Mood and Feelings Questionnaire depressive-symptom sum",
      "Age 23 standalone psychological distress outcome", "Age 23 standalone anxiety outcome", "Age 23 standalone depression outcome", "Age 23 standalone loneliness outcome", "Age 23 standalone wellbeing outcome",
      "Positive school wellbeing", "Positive schoolwork wellbeing", "Positive friend wellbeing", "Engagement/effort", "School engagement", "Reverse-scored disengagement",
      "Primary executive-control proxy used in the main models",
      "Alternative precision proxy sensitivity", "Alternative precision proxy sensitivity",
      "Parent conflict composite", "Sleep problem proxy", "Peer school trouble proxy", "Social support",
      "Exploratory genetic-score moderator; reverse-scored because lower raw values indicate higher depression liability in the CLS release.",
      "Exploratory genetic-score moderator for major depressive disorder liability.",
      rep("Exploratory genetic-score moderator", 10)
    )
  )
}

build_mcs_dataset <- function(config) {
  mcs_dir <- config$raw_data$mcs_dir
  pc_cols <- paste0("pc", 1:20)

  age11_derived <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs5_cm_derived.sav"),
    c("MCSID", "ECNUM00", "EEMOTION", "EPEER", "EHYPER", "EEBDTOT")
  ) %>%
    filter_mcs_first_cm("ECNUM00") %>%
    collapse_by_id("MCSID")

  age11_interview <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs5_cm_interview.sav"),
    c(
      "MCSID", "ECNUM00", "ECQ10A00", "ECQ10D00", "ECQ10E00", "ECQ29X00", "ECQ35X00",
      "ECQ36X00", "ECQ04X00", "ECQ17X00", "ECQ18X00", "ECQ48X00", "ECQ76X00"
    )
  ) %>%
    filter_mcs_first_cm("ECNUM00") %>%
    collapse_by_id("MCSID")

  age14_derived <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs6_cm_derived.sav"),
    c("MCSID", "FCNUM00", "FEMOTION", "FPEER", "FHYPER", "FCGTQOFDM", "FCGTDELAY", "FCGTRISKA", "FCGTRISKT", "FDCE0600")
  ) %>%
    filter_mcs_first_cm("FCNUM00") %>%
    collapse_by_id("MCSID")

  age14_interview <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs6_cm_interview.sav"),
    c(
      "MCSID", "FCNUM00", "FCCSEX00", "FCCAGE00", "FCSCHL00", "FCSCWK00", "FCFRNS00",
      "FCSCBE00", "FCSINT00", "FCSCWA00", "FCRJOY00", "FCMNWO00",
      paste0("FCMDS", LETTERS[1:13], "00"),
      "FCQUAM00", "FCQUAF00", "FCSLLN00", "FCPETR00", "FCSAFF00"
    )
  ) %>%
    filter_mcs_first_cm("FCNUM00") %>%
    collapse_by_id("MCSID")

  age14_family <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs6_family_derived.sav"),
    c("MCSID", "FOECDUK0")
  ) %>%
    collapse_by_id("MCSID")

  age14_cognitive <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs6_cm_cognitive_assessment.sav"),
    c("MCSID", "FCNUM00", "FCWRDSC")
  ) %>%
    filter_mcs_first_cm("FCNUM00") %>%
    collapse_by_id("MCSID")

  age17_derived <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs7_cm_derived.sav"),
    c("MCSID", "GCNUM00", "GDCKESSL", "GDWEMWBS", "GDCCONSC")
  ) %>%
    filter_mcs_first_cm("GCNUM00") %>%
    collapse_by_id("MCSID")

  age17_interview <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs7_cm_interview.sav"),
    c("MCSID", "GCNUM00", "GCWWOP00", "GCRJOY00", "GCSPFD00", "GCVICG00", "GCSQLT00", "GCCAPL00")
  ) %>%
    filter_mcs_first_cm("GCNUM00") %>%
    collapse_by_id("MCSID")

  age23_derived <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs8_23y_cm_derived.sav"),
    c("mcsid", "hcnum00", "hdageint", "hdcntry", "hdrgn", "hdimd", "hdkessler6", "hdwemwb", "hdloneliness", "hdgad2", "hdphq2")
  ) %>%
    filter_mcs_first_cm("hcnum00") %>%
    rename(MCSID = mcsid) %>%
    collapse_by_id("MCSID")

  prs <- read_mcs_sav_any(
    fs::path(mcs_dir, "cls_mcs_pgi_final_v1_1_mcsid.sav"),
    c(
      "mcsid", "pnum", pc_cols,
      "mcs_depression_p1_v1", "mcs_mdd_p1_v1", "mcs_anxiety_p1_v1", "mcs_adhd_p1_v1", "mcs_externalising_p5e08_v1",
      "mcs_cognition_p1_v1", "mcs_asd_p1_v1", "mcs_agreeableness_p1_v1", "mcs_consc_p1_v1",
      "mcs_extraversion_p1_v1", "mcs_neuroticism_p1_v1", "mcs_openness_p1_v1"
    )
  ) %>%
    filter_mcs_first_pnum("pnum") %>%
    rename(MCSID = mcsid) %>%
    collapse_by_id("MCSID") %>%
    mutate(across(all_of(pc_cols), clean_numeric))

  family_design <- read_mcs_sav_any(
    fs::path(mcs_dir, "mcs_longitudinal_family_file.sav"),
    c("MCSID", "PTTYPE2", "SPTN00", "FOVWT1", "FOVWT2", "GOVWT1", "GOVWT2")
  ) %>%
    collapse_by_id("MCSID") %>%
    mutate(
      mcs_design_stratum = na_if(clean_numeric(PTTYPE2), -1),
      mcs_sampling_stratum = mcs_pttype2_factor(mcs_design_stratum),
      mcs_design_psu = as.factor(clean_character(SPTN00)),
      mcs6_weight_uk = na_if(clean_numeric(FOVWT2), -1),
      mcs6_weight_country = na_if(clean_numeric(FOVWT1), -1),
      mcs7_weight_uk = na_if(clean_numeric(GOVWT2), -1),
      mcs7_weight_country = na_if(clean_numeric(GOVWT1), -1)
    ) %>%
    select(MCSID, mcs_design_stratum, mcs_sampling_stratum, mcs_design_psu, mcs6_weight_uk, mcs6_weight_country, mcs7_weight_uk, mcs7_weight_country)

  age11 <- age11_derived %>%
    left_join(age11_interview, by = "MCSID") %>%
    mutate(
      wave = "age11",
      age_years = 11,
      internalising_proxy = clean_numeric(EEMOTION) + clean_numeric(EPEER),
      distress_proxy = internalising_proxy,
      field_index = row_mean_standardised(
        pick(everything()),
        c("ECQ10A00", "ECQ10D00", "ECQ10E00", "ECQ29X00", "ECQ35X00", "ECQ36X00", "ECQ04X00", "ECQ48X00")
      ),
      precision_index = NA_real_,
      context_index = NA_real_
    )

  age14 <- age14_derived %>%
    left_join(age14_interview, by = "MCSID") %>%
    left_join(age14_cognitive, by = "MCSID") %>%
    left_join(age14_family, by = "MCSID") %>%
    mutate(
      wave = "age14",
      age_years = 14,
      sex = clean_numeric(FCCSEX00),
      age_baseline = clean_numeric(FCCAGE00),
      ethnicity_6cat = mcs_ethnicity_6cat(clean_numeric(FDCE0600)),
      internalising_proxy = clean_numeric(FEMOTION) + clean_numeric(FPEER),
      distress_proxy = internalising_proxy,
      school_happy_pos = -clean_numeric(FCSCHL00),
      schoolwork_happy_pos = -clean_numeric(FCSCWK00),
      friends_happy_pos = -clean_numeric(FCFRNS00),
      try_best_pos = -clean_numeric(FCSCBE00),
      school_interest_pos = -clean_numeric(FCSINT00),
      school_waste_reverse_pos = clean_numeric(FCSCWA00),
      read_enjoyment_pos = -clean_numeric(FCRJOY00),
      field_index = row_mean_standardised(
        pick(everything()),
        c(
          "school_happy_pos", "schoolwork_happy_pos", "friends_happy_pos",
          "try_best_pos", "school_interest_pos", "school_waste_reverse_pos"
        )
      ),
      quality_decision_pos = clean_numeric(FCGTQOFDM),
      delay_aversion_pos = -clean_numeric(FCGTDELAY),
      hyperactivity_reverse_pos = -clean_numeric(FHYPER),
      mind_on_schoolwork_reverse_pos = -clean_numeric(FCMNWO00),
      word_activity_pos = clean_numeric(FCWRDSC),
      precision_index = zscore(quality_decision_pos),
      precision_delay_index = zscore(delay_aversion_pos),
      precision_word_index = zscore(word_activity_pos),
      family_conflict = row_mean_standardised(pick(everything()), c("FCQUAM00", "FCQUAF00")),
      sleep_risk = zscore(FCSLLN00),
      peer_school_trouble = zscore(FCPETR00),
      social_support = zscore(FCSAFF00),
      income_quintile = na_if(clean_numeric(FOECDUK0), -1),
      income_quintile_z = zscore(income_quintile),
      mfq_age14_sum = row_sum_with_min(
        pick(everything()),
        paste0("FCMDS", LETTERS[1:13], "00"),
        min_non_missing = 10,
        transform = function(x) x - 1
      ),
      mfq_age14_z = zscore(mfq_age14_sum)
    )

  age17 <- age17_derived %>%
    left_join(age17_interview, by = "MCSID") %>%
    mutate(
      wave = "age17",
      age_years = 17,
      psychological_distress = clean_numeric(GDCKESSL),
      wellbeing = clean_numeric(GDWEMWBS),
      distress_proxy = psychological_distress,
      optimism_future_pos = clean_numeric(GCWWOP00),
      read_enjoyment_pos = -clean_numeric(GCRJOY00),
      friends_time_pos = -clean_numeric(GCSPFD00),
      future_alignment_pos = -clean_numeric(GCCAPL00),
      field_index = row_mean_standardised(
        pick(everything()),
        c("optimism_future_pos", "read_enjoyment_pos", "friends_time_pos", "future_alignment_pos")
      ),
      precision_index = zscore(GDCCONSC),
      context_index = row_mean_standardised(pick(everything()), c("GCVICG00", "GCSQLT00"))
    )

  age23 <- age23_derived %>%
    mutate(
      wave = "age23",
      age_years = 23,
      psychological_distress = clean_numeric(hdkessler6),
      anxiety_symptoms = clean_numeric(hdgad2),
      depressive_symptoms = clean_numeric(hdphq2),
      loneliness = clean_numeric(hdloneliness),
      wellbeing = clean_numeric(hdwemwb),
      distress_proxy = psychological_distress,
      field_index = NA_real_,
      precision_index = NA_real_,
      context_index = zscore(hdimd)
    )

  long <- bind_rows(
    age11 %>% select(MCSID, wave, age_years, distress_proxy, field_index, precision_index, context_index),
    age14 %>% select(MCSID, wave, age_years, distress_proxy, field_index, precision_index, family_conflict, sleep_risk, peer_school_trouble, social_support, income_quintile, mfq_age14_sum, sex, age_baseline),
    age17 %>% select(MCSID, wave, age_years, distress_proxy, field_index, precision_index, context_index),
    age23 %>% select(MCSID, wave, age_years, distress_proxy, field_index, precision_index, context_index)
  ) %>%
    group_by(wave) %>%
    mutate(distress_wave_z = zscore(distress_proxy)) %>%
    ungroup() %>%
    mutate(
      distress_pooled_z = zscore(distress_proxy),
      distress_z = distress_pooled_z
    )

  baseline <- age14 %>%
    select(
      MCSID, sex, age_baseline, ethnicity_6cat, distress_proxy, field_index, precision_index,
      precision_delay_index, precision_word_index,
      family_conflict, sleep_risk, peer_school_trouble, social_support,
      income_quintile, income_quintile_z, mfq_age14_sum, mfq_age14_z
    ) %>%
    mutate(
      baseline_distress_raw = distress_proxy,
      baseline_distress_z = zscore(distress_proxy)
    ) %>%
    left_join(family_design, by = "MCSID") %>%
    rename(sampling_stratum = mcs_sampling_stratum) %>%
    left_join(prs, by = "MCSID") %>%
    mutate(
      across(all_of(pc_cols), zscore),
      depression_prs = zscore(-mcs_depression_p1_v1),
      mdd_prs = zscore(mcs_mdd_p1_v1),
      anxiety_prs = zscore(mcs_anxiety_p1_v1),
      adhd_prs = zscore(mcs_adhd_p1_v1),
      externalising_prs = zscore(mcs_externalising_p5e08_v1),
      cognition_pgi = zscore(mcs_cognition_p1_v1),
      asd_pgi = zscore(mcs_asd_p1_v1),
      agreeableness_pgi = zscore(mcs_agreeableness_p1_v1),
      conscientiousness_pgi = zscore(mcs_consc_p1_v1),
      extraversion_pgi = zscore(mcs_extraversion_p1_v1),
      neuroticism_pgi = zscore(mcs_neuroticism_p1_v1),
      openness_pgi = zscore(mcs_openness_p1_v1)
    ) %>%
    select(
      -distress_proxy, -mcs_depression_p1_v1, -mcs_mdd_p1_v1, -mcs_anxiety_p1_v1,
      -mcs_adhd_p1_v1, -mcs_externalising_p5e08_v1,
      -mcs_cognition_p1_v1, -mcs_asd_p1_v1, -mcs_agreeableness_p1_v1,
      -mcs_consc_p1_v1, -mcs_extraversion_p1_v1, -mcs_neuroticism_p1_v1,
      -mcs_openness_p1_v1
    )

  predictive_age17 <- baseline %>%
    inner_join(
      age17 %>%
        transmute(
          MCSID,
          followup_wave = "age17",
          followup_age = 17,
          followup_distress_z = zscore(distress_proxy),
          followup_wellbeing_z = zscore(wellbeing)
        ),
      by = "MCSID"
    ) %>%
    mutate(
      exposure_wave = "age14",
      analysis_weight_uk = mcs6_weight_uk,
      analysis_weight_country = mcs6_weight_country,
      design_stratum = mcs_design_stratum,
      design_psu = mcs_design_psu
    )

  predictive_age23 <- baseline %>%
    inner_join(
      age23 %>%
        transmute(
          MCSID,
          followup_wave = "age23",
          followup_age = 23,
          followup_distress_z = zscore(distress_proxy),
          followup_anxiety_z = zscore(anxiety_symptoms),
          followup_depression_z = zscore(depressive_symptoms),
          followup_loneliness_z = zscore(loneliness),
          followup_wellbeing_z = zscore(wellbeing)
        ),
      by = "MCSID"
    ) %>%
    mutate(
      exposure_wave = "age14",
      analysis_weight_uk = mcs6_weight_uk,
      analysis_weight_country = mcs6_weight_country,
      design_stratum = mcs_design_stratum,
      design_psu = mcs_design_psu
    )

  predictive <- bind_rows(predictive_age17, predictive_age23)

  reliability_inputs <- list(
    field = age14 %>%
      select(
        MCSID,
        school_happy_pos,
        schoolwork_happy_pos,
        friends_happy_pos,
        try_best_pos,
        school_interest_pos,
        school_waste_reverse_pos
      ),
    precision = age14 %>%
      select(
        MCSID,
        quality_decision_pos
      )
  )

  wide <- predictive_age23 %>%
    left_join(
      age17 %>%
        transmute(
          MCSID,
          gdckessl_17 = psychological_distress,
          gdwemwbs_17 = wellbeing,
          age17_age = age_years
        ),
      by = "MCSID"
    ) %>%
    left_join(
      age23 %>%
        transmute(
          MCSID,
          hdkessler6_23 = psychological_distress,
          anxiety_23 = anxiety_symptoms,
          depression_23 = depressive_symptoms,
          loneliness_23 = loneliness,
          wellbeing_23 = wellbeing
        ),
      by = "MCSID"
    ) %>%
    mutate(
      age17_age = first_non_missing(age17_age, age_baseline + 3),
      GOVWT2 = mcs7_weight_uk,
      PTTYPE2 = mcs_design_stratum,
      SPTN00 = as.character(mcs_design_psu)
    )

  list(
    long = long,
    predictive = predictive,
    predictive_age17 = predictive_age17,
    predictive_age23 = predictive_age23,
    wide = wide,
    reliability_inputs = reliability_inputs,
    family_design = family_design,
    survey_weight_note = "MCS inferential models use the survey weight from the exposure wave together with PTTYPE2 and SPTN00. Age-14 exposure models use FOVWT2; age-17 exposure models, where used, should use GOVWT2."
  )
}
