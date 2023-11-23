# targets.R

# setup -------------------------------------------------------------------

#pipeline packages
suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
})

# target options
tar_option_set(
  packages = c(
    "tidyverse",
    "patchwork",
    "scales"
  ),
  format = "qs",
  controller = crew::crew_controller_local(workers = 4)
)

# source functions
tar_source()


# targets -----------------------------------------------------------------

list(

  # blots -------------------------------------------------------------------

  tar_target(
    blots_file,
    raw_data_path("blots.csv"),
    format = "file_fast"
  ),
  tar_target(
    blots_raw,
    readr::read_csv(blots_file, show_col_types = FALSE)
  ),
  tar_target(
    blots_clean,
    clean_blots(blots_raw)
  ),
  tar_target(
    blots_filtered,
    filter_blots(blots_clean)
  ),
  tar_target(
    blots_index,
    index_blots(blots_filtered)
  ),
  tar_map(
    values = tibble::tribble(
      ~experiment,   ~protein,     ~name,              ~paired, ~comp,
      "ar",          "α-SMA",      "ar_sma",           FALSE,    rlang::sym("comps_ar"),
      "bleo",        "MCT1",       "bleo_mct1",        FALSE,    rlang::sym("comps_bleo"),
      "bleo",        "MCT4",       "bleo_mct4",        FALSE,    rlang::sym("comps_bleo"),
      "bleo",        "α-SMA",      "bleo_sma",         FALSE,    rlang::sym("comps_bleo"),
      "comp",        "MCT1",       "comp_mct1",        FALSE,    NULL,
      "comp",        "MCT4",       "comp_mct4",        FALSE,    NULL,
      "comp",        "α-SMA",      "comp_sma",         FALSE,    NULL,
      "diclofenac",  "α-SMA",      "diclofenac_sma",   FALSE,    rlang::sym("comps_nsaid"),
      "dual",        "AKT",        "dual_akt",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "AMPK",       "dual_ampk",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "CNN1",       "dual_cnn1",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Col1a1",     "dual_col1a1",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Col3a1",     "dual_col3a1",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "ERK",        "dual_erk",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "FN1",        "dual_fn1",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Kla",        "dual_kla",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "MCT1",       "dual_mct1",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "MCT4",       "dual_mct4",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Smad3",      "dual_smad3",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pAKT",       "dual_pakt",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pAMPK",      "dual_pampk",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pERK",       "dual_perk",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pSmad3",     "dual_psmad3",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "α-SMA",      "dual_sma",         FALSE,    rlang::sym("comps_dual"),
      "ipf-lf",      "Col1a1",     "ipf_lf_col1a1",    FALSE,    rlang::sym("comps_azd_2"),
      "ipf-lf",      "α-SMA",      "ipf_lf_sma",       FALSE,    rlang::sym("comps_azd_2"),
      "ipf",         "MCT1",       "ipf_mct1",         FALSE,     NULL,
      "ipf",         "MCT4",       "ipf_mct4",         FALSE,     NULL,
      "ipf",         "α-SMA",      "ipf_sma",          FALSE,     NULL,
      "lactate",     "Col1a1",     "lactate_col1a1",   FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "FN1",        "lactate_fn1",      FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "MCT1",       "lactate_mct1",     FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "MCT4",       "lactate_mct4",     FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "α-SMA",      "lactate_sma",      FALSE,    rlang::sym("comps_lactate"),
      "sirna",       "MCT1",       "sirna_mct1",       FALSE,    rlang::sym("comps_sirna"),
      "sirna",       "MCT4",       "sirna_mct4",       FALSE,    rlang::sym("comps_sirna"),
      "sirna",       "α-SMA",      "sirna_sma",        FALSE,    rlang::sym("comps_sirna"),
      "six",         "HIF-1α",     "six_hif1a",        FALSE,    rlang::sym("comps_azd_2"),
      "untreated",   "Col1a1",     "untreated_col1a1", FALSE,    NULL,
      "untreated",   "FN1",        "untreated_fn1",    FALSE,    NULL,
      "untreated",   "MCT1",       "untreated_mct1",   FALSE,    NULL,
      "untreated",   "MCT4",       "untreated_mct4",   FALSE,    NULL,
      "untreated",   "α-SMA",      "untreated_sma",    FALSE,    NULL
    ),
    names = name,
    tar_target(
      blot_norm,
      norm_blot(blots_filtered, experiment, protein)
    ),
    tar_target(
      blot_stats,
      analyze_blot(blot_norm, paired, comp)
    ),
    tar_target(
      blot_plot,
      plot_blot(blot_norm, blot_stats, title = protein),
      format = "rds"
    ),
    # tar_target(
    #   blot_img,
    #   write_plot(blot_plot, name, path = "analysis/figures/blots"),
    #   format = "file"
    # ),
    NULL
  ),

  NULL
)
