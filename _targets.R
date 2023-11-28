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
      "ipf",         "MCT1",       "ipf_mct1",         FALSE,    NULL,
      "ipf",         "MCT4",       "ipf_mct4",         FALSE,    NULL,
      "ipf",         "α-SMA",      "ipf_sma",          FALSE,    NULL,
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
    NULL
  ),
  tar_target(
    blots_phospho,
    phospho_ratio(blots_filtered)
  ),
  tar_map(
    values = tibble::tribble(
      ~experiment,   ~protein,     ~name,              ~paired, ~comp,
      "dual",        "AKT",        "dual_p_akt",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "AMPK",       "dual_p_ampk",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "ERK",        "dual_p_erk",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Smad3",      "dual_p_smad3",     FALSE,    rlang::sym("comps_azd_2")
    ),
    names = name,
    tar_target(
      blot_norm,
      norm_blot(blots_phospho, experiment, protein)
    ),
    tar_target(
      blot_stats,
      analyze_blot(blot_norm, paired, comp)
    ),
    tar_target(
      blot_plot,
      plot_blot(blot_norm, blot_stats, title = paste0("p", protein, " / ", protein)),
      format = "rds"
    ),
    NULL
  ),

  # contraction -------------------------------------------------------------

  tar_target(
    contraction_file,
    raw_data_path("contraction.csv"),
    format = "file_fast"
  ),
  tar_target(
    contraction_raw,
    readr::read_csv(contraction_file, show_col_types = FALSE)
  ),
  tar_target(
    contraction_clean,
    clean_contraction(contraction_raw)
  ),
  tar_target(
    contraction_stats,
    onefactor(contraction_clean, "frac_ctl", comps = comps_azd_1)
  ),
  tar_target(
    contraction_plot,
    plot_two_factor(
      contraction_clean,
      contraction_stats,
      "treatment",
      "frac_ctl",
      ytitle = "Contracted area\n(TGFβ/Ctl)"
    ) +
      ggplot2::facet_wrap(
        ggplot2::vars("experiment"),
        strip.position = "top"
      ) +
      ggplot2::guides(fill = "none") +
      ggplot2::theme(
        strip.text = ggplot2::element_blank()
      ),
    format = "rds"
  ),

  # plates ------------------------------------------------------------------

  tar_map(
    values = list(
      x = c("picogreen", "lactate", "nad_protein", "nad_nad", "glucose"),
      cf = c(5129, 20000, 300, 900, 1110000),
      bywell = c(TRUE, TRUE, FALSE, FALSE, TRUE)
    ),
    names = x,
    tar_target(
      plates_file,
      raw_data_path(stringr::str_c(x, ".csv")),
      format = "file_fast"
    ),
    tar_target(
      plates_raw,
      readr::read_csv(plates_file, show_col_types = FALSE)
    ),
    tar_target(
      plates_clean,
      clean_plates(plates_raw)
    ),
    tar_target(
      plates_std_curves,
      make_std_curves(plates_clean),
      format = "rds"
    ),
    tar_target(
      plates_std_plots,
      write_plot_list(
        plates_std_curves$plots,
        plates_std_curves$title,
        stringr::str_c(x, "_std_curves")
      ),
      format = "file_fast"
    ),
    tar_target(
      plates_interp,
      interp_data(plates_clean, plates_std_curves)
    ),
    tar_target(
      plates_conc,
      clean_plates_interp(plates_interp, cf)
    ),
    NULL
  ),

  # extracellular -----------------------------------------------------------

  tar_map(
    values = list(
      df = rlang::syms(c("plates_conc_lactate", "plates_conc_glucose")),
      cf = rlang::syms(c("plates_conc_picogreen", "plates_conc_picogreen")),
      names = c("lactate", "glucose")
    ),
    names = names,
    tar_target(
      plates_norm,
      norm_plates(df, cf)
    )
  ),
  tar_map(
    values = list(
      names = c(
        "lactate_sirna",
        "lactate_azd",
        "lactate_ar",
        "lactate_nsaid",
        "glucose_azd"
      ),
      norm = rlang::syms(c(
        "plates_norm_lactate",
        "plates_norm_lactate",
        "plates_norm_lactate",
        "plates_norm_lactate",
        "plates_norm_glucose"
      )),
      exp = list(
        "sirna",
        "dual",
        c("dual", "ar"),
        "diclofenac",
        "dual"
      ),
      treat = list(
        NULL,
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AR", "VB", "AR/VB"),
        NULL,
        c("Veh", "AZD", "VB", "AZD/VB")
      ),
      comps = rlang::syms(c(
        "comps_sirna",
        "comps_azd_2",
        "comps_ar_2",
        "comps_nsaid",
        "comps_azd_2"
      )),
      ytitle = c(
        "Lactate (pmol/cell)",
        "Lactate (pmol/cell)",
        "Lactate (pmol/cell)",
        "Lactate (pmol/cell)",
        "Glucose (pmol/cell)"
      )
    ),
    names = names,
    tar_target(
      plate_filter,
      filter_plates(norm, exp = exp, treat = treat)
    ),
    tar_target(
      plate_stats,
      twofactor(plate_filter, "norm", comps = comps)
    ),
    tar_target(
      plate_plot,
      plot_two_factor(
        plate_filter,
        plate_stats,
        "treatment",
        "norm",
        ytitle = ytitle
      ),
      format = "rds"
    )
  ),

  # picogreen ---------------------------------------------------------------

  tar_map(
    values = list(
      names = c(
        "sirna",
        "azd",
        "ar",
        "nsaid"
      ),
      exp = list(
        "sirna",
        "dual",
        c("dual", "ar"),
        "diclofenac"
      ),
      treat = list(
        NULL,
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AR", "VB", "AR/VB"),
        NULL
      ),
      comps = rlang::syms(c(
        "comps_sirna",
        "comps_azd_2",
        "comps_ar_2",
        "comps_nsaid"
      ))
    ),
    names = names,
    tar_target(
      pg_filter,
      filter_plates(plates_conc_picogreen, exp = exp, treat = treat)
    ),
    tar_target(
      pg_stats,
      twofactor(pg_filter, "conc", comps = comps)
    ),
    tar_target(
      pg_plot,
      plot_two_factor(
        pg_filter,
        pg_stats,
        "treatment",
        "conc",
        ytitle = "Cell number"
      ) +
        ggplot2::facet_wrap(
          ggplot2::vars("experiment")
        ) +
        ggplot2::scale_y_continuous(
          labels = scales::label_number(scale_cut = scales::cut_short_scale()),
          breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
          expand = ggplot2::expansion(c(0, 0)),
          limits = nice_limits
        ) +
        ggplot2::theme(
          strip.text = ggplot2::element_blank()
        ),
      format = "rds"
    )
  ),

  # analysis ----------------------------------------------------------------

  tar_quarto(
    blots_analysis,
    path = "analysis/blots.qmd"
  ),
  tar_quarto(
    extracellular_analysis,
    path = "analysis/extracellular.qmd"
  ),

  NULL
)
