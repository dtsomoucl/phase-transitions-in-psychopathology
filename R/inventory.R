### DT --> Data inventory builders for the local ABCD and MCS assets and the selected analytic measures.
abcd_selected_sources <- function(config) {
  tibble(
    cohort = "ABCD",
    source_type = "analytic_source",
    label = c(
      "Parent demographics",
      "CBCL outcomes",
      "BIS/BAS field proxy",
      "School engagement",
      "NIH Toolbox cognition",
      "Executive deficits",
      "Family environment",
      "Sleep timing",
      "Cyberbullying",
      "Neighborhood ADI"
    ),
    path = fs::path(
      config$raw_data$abcd_dir,
      c(
        "ab_p_demo.tsv",
        "mh_p_cbcl.tsv",
        "mh_y_bisbas.tsv",
        "fc_y_srpf.tsv",
        "nc_y_nihtb.tsv",
        "nc_p_bdefs.tsv",
        "fc_y_fes.tsv",
        "ph_y_mctq.tsv",
        "mh_y_cb.tsv",
        "le_l_adi.tsv"
      )
    )
  )
}

mcs_selected_sources <- function(config) {
  tibble(
    cohort = "MCS",
    source_type = "analytic_source",
    label = c(
      "Age 11 derived",
      "Age 11 CM interview",
      "Age 14 derived",
      "Age 14 CM interview",
      "Age 14 family derived",
      "Longitudinal family file",
      "Age 17 derived",
      "Age 17 CM interview",
      "Age 23 derived",
      "PGI/PRS protected file",
      "Longitudinal dictionary workbook"
    ),
    path = c(
      fs::path(config$raw_data$mcs_dir, "mcs5_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs5_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs6_family_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs_longitudinal_family_file.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs7_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs7_cm_interview.sav"),
      fs::path(config$raw_data$mcs_dir, "mcs8_23y_cm_derived.sav"),
      fs::path(config$raw_data$mcs_dir, "cls_pgi_v1_mcs_protect.sav"),
      fs::path(config$raw_data$mcs_dir, "All_data/Sweep6/mrdoc/excel/mcs_longitudinal_data_dictionary.xlsx")
    )
  )
}

file_inventory <- function(base_dir, cohort) {
  fs::dir_ls(base_dir, recurse = TRUE, type = "file") %>%
    tibble(path = .) %>%
    mutate(
      cohort = cohort,
      extension = fs::path_ext(path),
      bytes = file.info(path)$size
    )
}

write_inventories <- function(config, abcd_mapping, mcs_mapping) {
  inventory_dir <- config$outputs$inventory_dir

  abcd_all <- file_inventory(config$raw_data$abcd_dir, "ABCD")
  mcs_all <- file_inventory(config$raw_data$mcs_dir, "MCS")

  selected_sources <- bind_rows(
    abcd_selected_sources(config),
    mcs_selected_sources(config)
  ) %>%
    mutate(
      exists = file.exists(path),
      bytes = file.info(path)$size
    )

  safe_write_csv(abcd_all, fs::path(inventory_dir, "abcd_file_inventory.csv"))
  safe_write_csv(mcs_all, fs::path(inventory_dir, "mcs_file_inventory.csv"))
  safe_write_csv(selected_sources, fs::path(inventory_dir, "selected_source_inventory.csv"))
  safe_write_csv(abcd_mapping, fs::path(inventory_dir, "abcd_variable_mapping.csv"))
  safe_write_csv(mcs_mapping, fs::path(inventory_dir, "mcs_variable_mapping.csv"))

  bind_rows(abcd_all, mcs_all) %>%
    count(cohort, extension, name = "n_files") %>%
    arrange(cohort, desc(n_files)) %>%
    safe_write_csv(fs::path(inventory_dir, "file_type_summary.csv"))

  invisible(
    list(
      abcd_all = abcd_all,
      mcs_all = mcs_all,
      selected_sources = selected_sources
    )
  )
}
