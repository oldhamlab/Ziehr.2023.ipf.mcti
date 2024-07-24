# targets.R

# setup -------------------------------------------------------------------

#pipeline packages
suppressPackageStartupMessages({
  library(targets)
  library(tarchetypes)
  library(patchwork)
})

# target options
tar_option_set(
  packages = c(
    "tidyverse",
    "patchwork",
    "scales",
    "SummarizedExperiment",
    "extrafont"
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
    blots_phospho,
    blot_ratios(blots_filtered)
  ),
  tar_target(
    blots_all,
    dplyr::bind_rows(blots_filtered, blots_phospho)
  ),
  tar_target(
    blots_index,
    index_blots(blots_all)
  ),
  tar_map(
    values = tibble::tribble(
      ~experiment,   ~protein,       ~name,              ~paired, ~comp,
      "ar",          "α-SMA",        "ar_sma",           FALSE,    rlang::sym("comps_ar"),
      "bleo",        "MCT1",         "bleo_mct1",        FALSE,    rlang::sym("comps_bleo"),
      "bleo",        "MCT4",         "bleo_mct4",        FALSE,    rlang::sym("comps_bleo"),
      "bleo",        "α-SMA",        "bleo_sma",         FALSE,    rlang::sym("comps_bleo"),
      "comp",        "MCT1",         "comp_mct1",        FALSE,    NULL,
      "comp",        "MCT4",         "comp_mct4",        FALSE,    NULL,
      "comp",        "α-SMA",        "comp_sma",         FALSE,    NULL,
      "diclofenac",  "α-SMA",        "diclofenac_sma",   FALSE,    rlang::sym("comps_nsaid"),
      "dual",        "AKT",          "dual_akt",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "AMPK",         "dual_ampk",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "CNN1",         "dual_cnn1",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Col1a1",       "dual_col1a1",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Col3a1",       "dual_col3a1",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "ERK",          "dual_erk",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "FN1",          "dual_fn1",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Kla",          "dual_kla",         FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "MCT1",         "dual_mct1",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "MCT4",         "dual_mct4",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "Smad3",        "dual_smad3",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pAKT",         "dual_pakt",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pAMPK",        "dual_pampk",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pERK",         "dual_perk",        FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pSmad3",       "dual_psmad3",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "α-SMA",        "dual_sma",         FALSE,    rlang::sym("comps_dual"),
      "ipf-lf",      "Col1a1",       "ipf_lf_col1a1",    FALSE,    rlang::sym("comps_azd_2"),
      "ipf-lf",      "α-SMA",        "ipf_lf_sma",       FALSE,    rlang::sym("comps_azd_2"),
      "ipf",         "MCT1",         "ipf_mct1",         FALSE,    NULL,
      "ipf",         "MCT4",         "ipf_mct4",         FALSE,    NULL,
      "ipf",         "α-SMA",        "ipf_sma",          FALSE,    NULL,
      "lactate",     "Col1a1",       "lactate_col1a1",   FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "FN1",          "lactate_fn1",      FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "MCT1",         "lactate_mct1",     FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "MCT4",         "lactate_mct4",     FALSE,    rlang::sym("comps_lactate"),
      "lactate",     "α-SMA",        "lactate_sma",      FALSE,    rlang::sym("comps_lactate"),
      "sirna",       "MCT1",         "sirna_mct1",       FALSE,    rlang::sym("comps_sirna"),
      "sirna",       "MCT4",         "sirna_mct4",       FALSE,    rlang::sym("comps_sirna"),
      "sirna",       "α-SMA",        "sirna_sma",        FALSE,    rlang::sym("comps_sirna"),
      "six",         "HIF-1α",       "six_hif1a",        FALSE,    rlang::sym("comps_azd_2"),
      "untreated",   "Col1a1",       "untreated_col1a1", FALSE,    NULL,
      "untreated",   "FN1",          "untreated_fn1",    FALSE,    NULL,
      "untreated",   "MCT1",         "untreated_mct1",   FALSE,    NULL,
      "untreated",   "MCT4",         "untreated_mct4",   FALSE,    NULL,
      "untreated",   "α-SMA",        "untreated_sma",    FALSE,    NULL,
      "dual",        "pAKT/AKT",     "dual_p_akt",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pAMPK/AMPK",   "dual_p_ampk",      FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pERK/ERK",     "dual_p_erk",       FALSE,    rlang::sym("comps_azd_2"),
      "dual",        "pSmad3/Smad3", "dual_p_smad3",     FALSE,    rlang::sym("comps_azd_2"),
      "lactyl",      "Kla",          "lactyl_kla",       FALSE,    rlang::sym("comps_azd_2"),
      "lactyl",      "H3",           "lactyl_h3",        FALSE,    rlang::sym("comps_azd_2"),
      "lactyl",      "Kla/H3",       "lactyl_kla_h3",    FALSE,    rlang::sym("comps_azd_2")
    ),
    names = name,
    tar_target(
      blot_norm,
      norm_blot(blots_all, experiment, protein)
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
  tar_map(
    values = list(
      experiment = list(
        "ipf",
        "bleo",
        "untreated",
        "sirna",
        "dual",
        c("dual", "ar"),
        "ipf-lf",
        "dual",
        "dual",
        "lactyl"
      ),
      protein = list(
        c("α-SMA", "MCT1", "MCT4"),
        c("α-SMA", "MCT1", "MCT4"),
        c("α-SMA", "MCT1", "MCT4"),
        c("α-SMA", "MCT1", "MCT4"),
        c("α-SMA", "Col1a1"),
        "α-SMA",
        c("α-SMA", "Col1a1"),
        c("pSmad3", "Smad3", "pSmad3/Smad3"),
        c("pERK", "ERK", "pERK/ERK"),
        c("Kla", "H3", "Kla/H3")
      ),
      treatment = list(
        NULL,
        "Veh",
        NULL,
        NULL,
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AR", "VB", "AR/VB"),
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AZD", "VB", "AZD/VB"),
        c("Veh", "AZD", "VB", "AZD/VB")
      ),
      condition = list(
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL,
        NULL
      ),
      comp = list(
        NULL,
        NULL,
        NULL,
        rlang::sym("comps_sirna"),
        rlang::sym("comps_azd_2"),
        rlang::sym("comps_ar_2"),
        rlang::sym("comps_azd_2"),
        rlang::sym("comps_azd_2"),
        rlang::sym("comps_azd_2"),
        rlang::sym("comps_azd_2")
      ),
      fn = list(
        group_ratio_ttest,
        group_ratio_ttest,
        group_ratio_ttest,
        group_twofactor,
        group_twofactor,
        group_twofactor,
        group_twofactor,
        group_twofactor,
        group_twofactor,
        group_twofactor
      ),
      nrow = list(1, 1, 1, NULL, NULL, NULL, NULL, NULL, NULL, NULL),
      ncol = list(NULL, NULL, NULL, 1, 1, 1, 1, 1, 1, 1),
      strip_pos = list(
        "top",
        "top",
        "top",
        "right",
        "right",
        "top",
        "right",
        "right",
        "right",
        "right"
      ),
      name = list(
        "fig_ipf",
        "fig_bleo",
        "fig_tgfb",
        "fig_sirna",
        "fig_azd",
        "fig_ar",
        "fig_ipf_lf",
        "fig_psmad3",
        "fig_perk",
        "fig_kla_h3"
      )
    ),
    names = name,
    tar_target(
      blot_norm,
      norm_blot(blots_all, experiment, protein, treatment, condition)
    ),
    tar_target(
      blot_stats,
      fn(blot_norm, "norm", comp)
    ),
    tar_target(
      blot_plot,
      plot_blot(blot_norm, blot_stats) +
        ggplot2::facet_wrap(
          ggplot2::vars(.data[["protein"]]),
          scales = "free_y",
          nrow = nrow,
          ncol = ncol,
          strip.position = strip_pos
        ),
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
    values = tibble::tribble(
      ~x, ~cf,
      "picogreen",    5129,
      "lactate",      20000,
      "nad_protein",  300,
      "nad_nad",      900,
      "glucose",      1110000,
      "nadp_dna",     5129,
      "nadp_nadp",    600
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
      filter_plates(plates_conc_picogreen, exp = exp, treat = treat) |>
        normalize("conc", "experiment")
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
        "conc_corr",
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

  # nad ---------------------------------------------------------------------

  tar_target(
    nad_clean,
    clean_nad(plates_conc_nad_protein, plates_conc_nad_nad)
  ),
  tar_target(
    nadp_clean,
    clean_nad(plates_conc_nadp_dna, plates_conc_nadp_nadp)
  ),
  tar_target(
    nad_nadp_clean,
    dplyr::bind_rows(nad_clean, nadp_clean)
  ),
  tar_map(
    values = list(
      names = c(
        "nad",
        "nad_norm",
        "nadh",
        "nadh_norm",
        "nadh_ratio",
        "nadp",
        "nadp_norm",
        "nadph",
        "nadph_norm",
        "nadph_ratio"
      ),
      measure = c(
        "NAD",
        "NAD_norm",
        "NADH",
        "NADH_norm",
        "NADH/NAD",
        "NADP",
        "NADP_norm",
        "NADPH",
        "NADPH_norm",
        "NADPH/NADP"
      ),
      ytitle = c(
        "NAD (pmol)",
        "NAD (pmol / μg protein)",
        "NADH (pmol)",
        "NADH (pmol / μg protein)",
        "NADH / NAD<sup>+</sup>",
        "NADP (pmol)",
        "NADP (fmol / cell)",
        "NADPH (pmol)",
        "NADPH (fmol / cell)",
        "NADPH / NADP<sup>+</sup>"
      )
    ),
    names = names,
    tar_target(
      nad_filter,
      dplyr::filter(nad_nadp_clean, measurement == measure)
    ),
    tar_target(
      nad_stats,
      twofactor(nad_filter, y = "value", comps = comps_azd_2, log = TRUE)
    ),
    tar_target(
      nad_norm,
      norm_nad(nad_filter)
    ),
    tar_target(
      nad_plot,
      plot_two_factor(
        nad_norm,
        nad_stats,
        "treatment",
        "value_corr",
        ytitle = ytitle
      ) +
        ggplot2::theme(axis.title.y.left = ggtext::element_markdown()),
      format = "rds"
    )
  ),

  # ros ---------------------------------------------------------------------

  tar_map(
    values = list(
      names = c("cellrox", "mitosox")
    ),
    names = names,
    tar_target(
      ros_file,
      raw_data_path(stringr::str_c(names, ".csv")),
      format = "file_fast"
    ),
    tar_target(
      ros_raw,
      readr::read_csv(ros_file, show_col_types = FALSE)
    ),
    tar_target(
      ros_clean,
      clean_ros(ros_raw)
    )
  ),
  tar_target(
    ros_clean,
    dplyr::bind_rows(ros_clean_cellrox, ros_clean_mitosox)
  ),
  tar_map(
    values = list(
      names = c("cellrox", "mitosox", "mitotracker", "ratio"),
      measure = c("CellROX", "MitoSOX", "MitoTracker", "Ratio"),
      ytitle = c(
        "CellROX",
        "MitoSOX",
        "MitoTracker",
        "MitoSOX (normalized)"
      )
    ),
    names = names,
    tar_target(
      ros_norm,
      norm_ros(ros_clean, measure)
    ),
    tar_target(
      ros_stats,
      seahorse_group_onefactor(ros_norm) |>
        dplyr::mutate(measurement = measure) |>
        dplyr::rename(group = "x")
    ),
    tar_target(
      ros_plot,
      plot_mice(
        ros_norm,
        ros_stats,
        measure = measure,
        x = "group",
        y = "value_corr",
        ytitle = ytitle
      ),
      format = "rds"
    ),
    NULL
  ),

  # bcecf -------------------------------------------------------------------

  tar_target(
    bcecf_file,
    raw_data_path("bcecf.csv"),
    format = "file_fast"
  ),
  tar_target(
    bcecf_data,
    readr::read_csv(bcecf_file, show_col_types = FALSE)
  ),
  tar_target(
    bcecf_raw,
    format_bcecf(bcecf_data)
  ),
  tar_target(
    bcecf_clean,
    clean_bcecf(bcecf_raw)
  ),
  tar_target(
    bcecf_stats,
    twofactor(bcecf_clean, "ratio", comps = comps_azd_2)
  ),
  tar_target(
    bcecf_plot,
    plot_bcecf(bcecf_clean, bcecf_stats)
  ),

  # seahorse ----------------------------------------------------------------

  tar_target(
    seahorse_mcti_meta_file,
    raw_data_path("seahorse_mcti.csv"),
    format = "file_fast"
  ),
  tar_target(
    seahorse_mcti_wells,
    format_wells(seahorse_mcti_meta_file)
  ),
  tar_target(
    seahorse_mcti_count_files,
    raw_data_path("seahorse_\\d{4}-\\d{2}-\\d{2}_mcti_counts\\.csv"),
    format = "file"
  ),
  tar_target(
    seahorse_mcti_counts,
    purrr::map(
      seahorse_mcti_count_files,
      \(x) seahorse::format_cells(x) |>
        dplyr::mutate(value = value / 1000)
    )
  ),
  tar_target(
    seahorse_mcti_data_files,
    raw_data_path("seahorse_\\d{4}-\\d{2}-\\d{2}_mcti\\.xlsx"),
    format = "file"
  ),
  tar_target(
    seahorse_mcti_raw,
    purrr::map2(
      seahorse_mcti_data_files,
      seahorse_mcti_counts,
      \(x, y) seahorse::Seahorse(
        path = x,
        wells = seahorse_mcti_wells,
        stages = seahorse::stages_mst,
        cells = y,
        bf = 2.4,
        cf = 0.410
      )
    )
  ),
  tar_target(
    seahorse_mcti_herd,
    seahorse::Herd(seahorse_mcti_raw)
  ),
  tar_target(
    seahorse_mcti_timeline_plot,
    plot_seahorse_timelines(seahorse_mcti_herd),
    format = "rds"
  ),
  tar_target(
    seahorse_mcti_summary_data,
    summarize_seahorse(seahorse_mcti_herd)
  ),
  tar_target(
    seahorse_mcti_summary_stats,
    seahorse_group_onefactor(seahorse_mcti_summary_data)
  ),
  tar_target(
    seahorse_mcti_summary_plot,
    plot_seahorse_summary(
      seahorse_mcti_summary_data,
      seahorse_mcti_summary_stats
    ),
    format = "rds"
  ),
  tar_target(
    seahorse_mcti_atp_data,
    summarize_seahorse_atp(seahorse_mcti_herd)
  ),
  tar_target(
    seahorse_mcti_atp_stats,
    seahorse_group_onefactor(seahorse_mcti_atp_data)
  ),
  tar_target(
    seahorse_mcti_atp_bars,
    plot_seahorse_atp_bars(
      seahorse_mcti_atp_data,
      seahorse_mcti_atp_stats
    ),
    format = "rds"
  ),
  tar_target(
    seahorse_mcti_atp_pheno,
    plot_seahorse_pheno(seahorse_mcti_atp_data),
    format = "rds"
  ),
  tar_target(
    seahorse_mcti_stress_data,
    summarize_seahorse_stress(seahorse_mcti_herd)
  ),
  tar_target(
    seahorse_mcti_stress_stats,
    seahorse_group_onefactor(seahorse_mcti_stress_data)
  ),
  tar_target(
    seahorse_mcti_stress_plot,
    plot_seahorse_stress(
      seahorse_mcti_stress_data,
      seahorse_mcti_stress_stats
    ),
    format = "rds"
  ),
  tar_target(
    seahorse_mcti_time_data_files,
    raw_data_path("mcti_timecourse\\.xlsx"),
    format = "file"
  ),
  tar_target(
    seahorse_mcti_time_raw,
    seahorse::Seahorse(
      path = seahorse_mcti_time_data_files,
      wells = seahorse_mcti_wells,
      bf = 2.4,
      cf = 0.410
    )
  ),
  tar_target(
    seahorse_mcti_time_timeline_plot,
    seahorse::plot(seahorse_mcti_time_raw, "rates"),
    format = "rds"
  ),

  # metab -------------------------------------------------------------------

  tar_target(
    metabolite_pathways_file,
    raw_data_path("metabolites.tab"),
    format = "file"
  ),
  tar_target(
    metabolite_pathways,
    read_pathways(metabolite_pathways_file)
  ),

  # metab mcti --------------------------------------------------------------

  tar_map(
    values = list(
      type = c("intra", "extra"),
      path = c(
        "myofib_azd3965-vb124_intra_2022-02-10",
        "myofib_azd3965-vb124_extra_2024-01-19"
      ),
      outs = list(5, 5)
    ),
    names = type,
    tar_target(
      metab_mcti_tar_file,
      raw_data_path(stringr::str_c(path, ".xlsx")),
      format = "file"
    ),
    tar_target(
      metab_mcti_tar_raw,
      format_metab_targeted(metab_mcti_tar_file)
    ),
    tar_target(
      metab_mcti_samples_file,
      raw_data_path(stringr::str_c(path, ".csv")),
      format = "file"
    ),
    tar_target(
      metab_mcti_samples,
      make_mcti_samples(metab_mcti_samples_file)
    ),
    tar_target(
      metab_mcti_tar_se,
      format_metab_mcti_tar(metab_mcti_tar_raw, metab_mcti_samples)
    ),
    tar_target(
      metab_mcti_tar_clean,
      metab_mcti_tar_se |>
        remove_missing_metab() |>
        correct_drift() |>
        quality_control() |>
        impute_missing() |>
        pqn() |>
        annot_metabs()
    ),
    tar_target(
      metab_mcti_tar_pca,
      plot_pca_metab(metab_mcti_tar_clean),
      format = "rds"
    ),
    tar_target(
      metab_mcti_tar_out,
      metab_mcti_tar_clean[, metab_mcti_tar_clean$replicate %nin% outs]
    ),
    tar_target(
      metab_mcti_tar_treated,
      metab_mcti_tar_out[, metab_mcti_tar_out$treatment %nin% c("None", "Med")]
    ),
    tar_target(
      metab_mcti_tar_pca_filter,
      plot_pca_metab(metab_mcti_tar_treated, show_reps = FALSE),
      format = "rds"
    ),
    tar_target(
      metab_mcti_tar_limma,
      fit_limma_metab(
        metab_mcti_tar_out,
        tgfb_main = (TGFβ.None + TGFβ.Veh + TGFβ.AZD + TGFβ.VB + TGFβ.AZD.VB) / 5 -
          (Ctl.None + Ctl.Veh + Ctl.AZD + Ctl.VB + Ctl.AZD.VB) / 5,
        azd_main = (TGFβ.AZD + Ctl.AZD) / 2 - (TGFβ.Veh + Ctl.Veh) / 2,
        vb_main = (TGFβ.VB + Ctl.VB) / 2 - (TGFβ.Veh + Ctl.Veh) / 2,
        dual_main = (TGFβ.AZD.VB + Ctl.AZD.VB) / 2 - (TGFβ.Veh + Ctl.Veh) / 2
      )
    ),
    tar_map(
      values = list(
        comparison = c("tgfb_main", "azd_main", "vb_main", "dual_main"),
        fills = list(
          c("TGFβ", "Ctl"),
          c("AZD", "Veh"),
          c("VB", "Veh"),
          c("AZD/VB", "Veh")
        ),
        title = c("TGFβ", "AZD", "VB", "AZD/VB"),
        filename = stringr::str_c("msea_", c("tgfb", "azd", "vb", "dual")),
        names = c("tgfb", "azd", "vb", "dual")
      ),
      names = names,
      tar_target(
        metab_mcti_tar_tt,
        top_table_metab(metab_mcti_tar_out, metab_mcti_tar_limma, comparison)
      ),
      tar_target(
        metab_mcti_tar_vol,
        plot_volcano(
          metab_mcti_tar_tt,
          x = "logFC",
          y = "adj.P.Val",
          stat = "t",
          label = "metabolite",
          fills = fills,
          title = title
        ),
        format = "rds"
      ),
      tar_target(
        metab_mcti_tar_msea,
        run_msea(metab_mcti_tar_tt, metabolite_pathways)
      ),
      tar_target(
        metab_mcti_tar_msea_table,
        plot_msea_table(metab_mcti_tar_msea, title = title, fills = fills),
        format = "rds"
      ),
      tar_target(
        msea_table_file,
        write_table(
          metab_mcti_tar_msea_table,
          path = "analysis/figures/msea/",
          stringr::str_c(filename, "_", type)
        ),
        format = "file"
      ),
      tar_target(
        msea_table_img,
        plot_image(msea_table_file),
        format = "rds"
      ),
      NULL
    ),
    tar_target(
      metab_mcti_tar_tt,
      dplyr::bind_rows(
        AZD = metab_mcti_tar_tt_azd,
        VB = metab_mcti_tar_tt_vb,
        `AZD/VB` = metab_mcti_tar_tt_dual,
        .id = "treatment"
      )
    ),
    tar_map(
      values = list(
        names = c("lactate", "proline")
      ),
      names = names,
      tar_target(
        metab_mcti_tar,
        plot_mois(metab_mcti_tar_treated, metab_mcti_tar_tt, names),
        format = "rds"
      )
    ),
    NULL
  ),

  # metab bleo --------------------------------------------------------------

  tar_target(
    metab_bleo_lung_tar_file,
    raw_data_path("bleo_lung_2022-11-22.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_bleo_plasma_tar_file,
    raw_data_path("bleo_plasma_2022-11-22.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_bleo_tar_samples,
    make_bleo_samples(mice_raw)
  ),
  tar_map(
    values = list(
      source = c("lung", "plasma"),
      input = rlang::syms(c(
        "metab_bleo_lung_tar_file",
        "metab_bleo_plasma_tar_file"
      ))
    ),
    names = source,
    tar_target(
      metab_bleo_tar_raw,
      format_metab_targeted(input)
    ),
    tar_target(
      metab_bleo_tar_se,
      format_metab_bleo_tar(metab_bleo_tar_raw, metab_bleo_tar_samples)
    ),
    tar_target(
      metab_bleo_tar_clean,
      metab_bleo_tar_se |>
        remove_missing_metab() |>
        correct_drift() |>
        quality_control() |>
        impute_missing() |>
        pqn() |>
        annot_metabs()
    ),
    tar_target(
      metab_bleo_tar_pca,
      plot_pca_metab(metab_bleo_tar_clean, batch = FALSE, show_reps = FALSE),
      format = "rds"
    ),
    tar_target(
      metab_bleo_tar_limma,
      fit_limma_metab(
        metab_bleo_tar_clean,
        batch = FALSE,
        bleo = Bleo - Ctl,
        azd = AZD - Bleo,
        vb = VB - Bleo
      )
    ),
    tar_map(
      values = list(
        comparison = c("bleo", "azd", "vb"),
        fills = list(
          c("Bleo", "Ctl"),
          c("AZD", "Bleo"),
          c("VB", "Bleo")
        ),
        title = c("Bleo", "AZD", "VB"),
        filename = stringr::str_c("msea_", c("bleo", "bleo_azd", "bleo_vb"))
      ),
      names = comparison,
      tar_target(
        metab_bleo_tar_tt,
        top_table_metab(
          metab_bleo_tar_clean,
          metab_bleo_tar_limma,
          comparison
        )
      ),
      tar_target(
        metab_bleo_tar_vol,
        plot_volcano(
          metab_bleo_tar_tt,
          x = "logFC",
          y = "adj.P.Val",
          stat = "t",
          label = "metabolite",
          fills = fills,
          title = title
        ),
        format = "rds"
      ),
      tar_target(
        metab_bleo_tar_msea,
        run_msea(metab_bleo_tar_tt, metabolite_pathways)
      ),
      tar_target(
        metab_bleo_tar_msea_table,
        plot_msea_table(metab_bleo_tar_msea, title = title, fills = fills),
        format = "rds"
      ),
      tar_target(
        msea_table_file,
        write_table(
          metab_bleo_tar_msea_table,
          path = "analysis/figures/msea/",
          stringr::str_c(filename, "_", source)
        ),
        format = "file"
      ),
      tar_target(
        msea_table_img,
        plot_image(msea_table_file),
        format = "rds"
      ),
      NULL
    ),
    NULL
  ),
  tar_target(
    metab_bleo_lacate_plot,
    plot_bleo_lactate(
      list(
        plasma = metab_bleo_tar_clean_plasma,
        lung = metab_bleo_tar_clean_lung
      )
    ),
    format = "rds"
  ),

  # mid ---------------------------------------------------------------------

  tar_target(
    qbias_file,
    raw_data_path("qbias-correction-factors"),
    format = "file"
  ),
  tar_target(
    qbias_ratios,
    format_qbias(qbias_file)
  ),
  tar_target(
    predicted_ratios,
    calculate_predicted_ratios(isotope_library)
  ),
  tar_target(
    correction_factors,
    calculate_correction_factors(qbias_ratios, predicted_ratios)
  ),
  tar_target(
    correction_matrices,
    make_correction_matrices(isotope_library)
  ),

  # mid mcti ----------------------------------------------------------------

  tar_target(
    mid_mcti_files,
    raw_data_path("\\d{4}-\\d{2}-\\d{2}_mcti_(glc2|gln5|lac3|glc6)\\.xlsx"),
    format = "file"
  ),
  tar_target(
    mid_mcti_mids,
    format_mid_mcti(mid_mcti_files, correction_factors)
  ),
  tar_target(
    mid_mcti_mids_corr,
    correct_mids(mid_mcti_mids, correction_matrices)
  ),
  tar_target(
    mid_mcti_stacked_data,
    format_mids_all(mid_mcti_mids_corr)
  ),
  tar_map(
    values = list(
      label = c("glc2", "gln5", "lac3", "glc6")
    ),
    names = label,
    tar_target(
      mid_mcti_plot,
      plot_mids_all(mid_mcti_stacked_data, label = label),
      format = "rds"
    ),
    tar_target(
      mid_mcti_img,
      write_plot(mid_mcti_plot, label, path = "analysis/figures/mids"),
      format = "file"
    ),
    NULL
  ),
  tar_map(
    values = list(
      label = c("gln5", "lac3", "glc6"),
      title = list(
        expression(paste("[U-"^13, "C"[5], "]-glutamine")),
        expression(paste("[U-"^13, "C"[3], "]-lactate")),
        expression(paste("[U-"^13, "C"[6], "]-glucose"))
      )
    ),
    names = label,
    tar_target(
      mid_mcti_moi_data,
      select_mid_mois(mid_mcti_mids_corr, label)
    ),
    tar_target(
      mid_mcti_moi_stats,
      group_twofactor(mid_mcti_moi_data, "labeled", comps_azd_2)
    ),
    tar_target(
      mid_mcti_moi_plot,
      plot_two_factor(
        mid_mcti_moi_data,
        mid_mcti_moi_stats,
        x = "treatment",
        y = "labeled",
        ytitle = "Labeled fraction",
        wrap = "metabolite",
        title = title
      ) +
        ggplot2::scale_y_continuous(
          breaks = seq(0, 1, 0.2),
          limits = c(0, 1),
          expand = ggplot2::expansion(c(0, 0))
        ),
      format = "rds"
    ),
    NULL
  ),
  tar_map(
    values = list(
      names = c(
        "glc_cit_m2_pyr",
        "glc_cit_m3_pyr",
        "glc_cit_m4_pyr",
        "glc_lac_pyr"
      ),
      top = alist(
        list(tracer = "glc6", metabolite = "CIT", isotope = "M2"),
        list(tracer = "glc6", metabolite = "CIT", isotope = "M3"),
        list(tracer = "glc6", metabolite = "CIT", isotope = "M4"),
        list(tracer = "glc6", metabolite = "LAC", isotope = "M3")
      ),
      bottom = alist(
        list(tracer = "glc6", metabolite = "PYR", isotope = "M3"),
        list(tracer = "glc6", metabolite = "PYR", isotope = "M3"),
        list(tracer = "glc6", metabolite = "PYR", isotope = "M3"),
        list(tracer = "glc6", metabolite = "PYR", isotope = "M3")
      ),
      ytitle = c(
        "CIT M2 / PYR M3",
        "CIT M3 / PYR M3",
        "CIT M4 / PYR M3",
        "LAC M3 / PYR M3"
      )
    ),
    names = names,
    tar_target(
      mid_ratios,
      calc_mid_ratio(mid_mcti_mids_corr, top, bottom)
    ),
    tar_target(
      mid_ratios_stats,
      twofactor(mid_ratios, "value", comps = comps_azd_2)
    ),
    tar_target(
      mid_ratios_plots,
      plot_two_factor(
        mid_ratios,
        mid_ratios_stats,
        x = "treatment",
        y = "value",
        ytitle = ytitle
      ) +
        ggplot2::scale_y_continuous(
          # breaks = seq(0, 1, 0.2),
          # limits = c(0, 1),
          expand = ggplot2::expansion(c(0, 0))
        ),
      format = "rds"
    ),
    NULL
  ),
  tar_map(
    values = list(
      names = c("gln5_proline"),
      tracer = c("gln5"),
      metabolite = c("PRO"),
      isotope = c("M5"),
      ytitle = c("M5 Mole fraction"),
      title = c("[U-^13^C~5~]-GLN → PRO")
    ),
    names = names,
    tar_target(
      mid_filter,
      filter_mids(
        mid_mcti_mids_corr,
        tracer,
        metabolite,
        isotope
      )
    ),
    tar_target(
      mid_stats,
      twofactor(
        mid_filter,
        "mid_corr",
        mixed = FALSE,
        comps = comps_azd_2
      )
    ),
    tar_target(
      mid_plot,
      plot_two_factor(
        mid_filter,
        mid_stats,
        x = "treatment",
        y = "mid_corr",
        ytitle = ytitle,
        title = title
      ) +
        ggplot2::theme(
          plot.title = ggtext::element_markdown(
            margin = ggplot2::margin(b = -5)
          )
        ),
      format = "rds"
    ),
    NULL
  ),
  tar_target(
    mid_tgf_plot,
    plot_tgf_mids(mid_mcti_mids_corr) +
      theme_patchwork(
        widths = ggplot2::unit(5.5, "in"),
        heights = ggplot2::unit(9, "in"),
        tags = NULL
      ),
    format = "rds"
  ),
  tar_target(
    fig05s,
    write_figures(mid_tgf_plot, "SF_05")
  ),

  # mid bleo ----------------------------------------------------------------

  tar_target(
    mid_bleo_files,
    raw_data_path("bleo_.*_isotope_2022-11-22.xlsx"),
    format = "file"
  ),
  tar_target(
    mid_bleo_mids,
    format_bleo_mids(mid_bleo_files, mice_raw)
  ),
  tar_target(
    mid_bleo_mids_corr,
    correct_mids(mid_bleo_mids, correction_matrices)
  ),
  tar_target(
    mids_bleo_stacked_data,
    format_mids_bleo_all(mid_bleo_mids_corr)
  ),
  tar_target(
    mid_bleo_plasma_glucose_plot,
    plot_plasma_glucose(mid_bleo_mids_corr),
    format = "rds"
  ),
  tar_target(
    mid_bleo_plasma_ratio,
    plot_plasma_ratio(mid_bleo_mids_corr),
    format = "rds"
  ),
  tar_map(
    values = list(
      source = c("plasma", "lung")
    ),
    names = source,
    tar_target(
      mids_bleo_plot,
      plot_mids_bleo_all(mids_bleo_stacked_data, source),
      format = "rds"
    ),
    tar_target(
      mid_bleo_img,
      write_plot(mids_bleo_plot, source, path = "analysis/figures/mids"),
      format = "file_fast"
    )
  ),
  tar_map(
    values = list(
      names = c(
        "lac3_glc6",
        "pg_glc6",
        "cit2_glc6",
        "lac3_pg",
        "cit2_pg"
      ),
      top = alist(
        list(tissue = "lung", metabolite = "LAC", isotope = "M3"),
        list(tissue = "lung", metabolite = "3PG", isotope = "M3"),
        list(tissue = "lung", metabolite = "CIT", isotope = "M2"),
        list(tissue = "lung", metabolite = "LAC", isotope = "M3"),
        list(tissue = "lung", metabolite = "CIT", isotope = "M2")
      ),
      bottom = alist(
        list(tissue = "plasma", metabolite = "GLC", isotope = "M6"),
        list(tissue = "plasma", metabolite = "GLC", isotope = "M6"),
        list(tissue = "plasma", metabolite = "GLC", isotope = "M6"),
        list(tissue = "lung", metabolite = "3PG", isotope = "M3"),
        list(tissue = "lung", metabolite = "3PG", isotope = "M3")
      ),
      ytitle = c(
        "LAC M3 / GLC M6",
        "3PG M3 / GLC M6",
        "CIT M2 / GLC M6",
        "LAC M3 / 3PG M3",
        "CIT M2 / 3PG M3"
      )
    ),
    names = names,
    tar_target(
      mid_ratios,
      calc_mid_ratio(mid_bleo_mids_corr, top, bottom) |>
        dplyr::group_by(group) |>
        wmo::remove_nested_outliers("value", remove = TRUE) |>
        identity()
    ),
    tar_target(
      mid_ratios_stats,
      stats_histo(mid_ratios, "ratio") |>
        dplyr::mutate(x = factor(group, levels = c("Ctl", "Bleo", "AZD", "VB")))
    ),
    tar_target(
      mid_ratios_plots,
      plot_one_factor(
        mid_ratios,
        mid_ratios_stats,
        x = "group",
        y = "value",
        ytitle = ytitle
      ),
      format = "rds"
    ),
    NULL
  ),

  # rnaseq ------------------------------------------------------------------

  tar_target(
    dds,
    count_rnaseq(rnaseq.lf.tgfb.mcti::se),
    format = "rds"
  ),
  tar_target(
    rnaseq_pca,
    vst_rnaseq(dds)
  ),
  tar_target(
    rnaseq_pca_plot,
    plot_rnaseq_pca(rnaseq_pca),
    format = "rds"
  ),
  tar_target(
    gsea_pathways,
    get_msigdb_pathways(category = "H")
  ),
  tar_target(
    tfea,
    run_tfea(dds) |>
      summarize_tf()
  ),
  tar_target(
    tfea_fit,
    fit_tfea(tfea)
  ),
  tar_map(
    values = list(
      names = c(
        "tgfb",
        "azd",
        "vb",
        "dual"
      ),
      con = list(
        c("group", "TGFβ.Veh", "Ctl.Veh"),
        c("group", "TGFβ.AZD", "TGFβ.Veh"),
        c("group", "TGFβ.VB", "TGFβ.Veh"),
        c("group", "TGFβ.AZD.VB", "TGFβ.Veh")
      ),
      fills = list(
        c("TGFβ", "Ctl"),
        c("AZD", "TGFβ"),
        c("VB", "TGFβ"),
        c("AZD/VB", "TGFβ")
      ),
      title = c(
        "TGFβ v. Ctl",
        "AZD v. Veh in TGFβ",
        "VB v. Veh in TGFβ",
        "AZD/VB v. Veh in TGFβ"
      )
    ),
    names = names,
    tar_target(
      deg,
      rnaseq_results(dds, con = con)
    ),
    tar_target(
      deg_vol,
      plot_vol_rnaseq(
        deg,
        x = "log2FoldChange",
        y = "padj",
        stat = "stat",
        label = "symbol",
        nudge_factor = 0.5,
        fills = fills,
        title = title
      ),
      format = "rds"
    ),
    tar_target(
      gsea,
      run_gsea(deg, gsea_pathways)
    ),
    tar_target(
      gsea_table,
      plot_msea_table(gsea, src = "HALLMARK", title = title, fills = fills),
      format = "rds"
    ),
    tar_target(
      gsea_table_file,
      write_table(
        gsea_table,
        path = "analysis/figures/gsea/",
        stringr::str_c("gsea_", names)
      ),
      format = "file"
    ),
    tar_target(
      gsea_table_img,
      plot_image(gsea_table_file),
      format = "rds"
    ),
    tar_target(
      tfea,
      index_tfea(tfea_fit, names)
    ),
    tar_target(
      tfea_vol,
      plot_volcano(
        tfea,
        x = "logFC",
        y = "adj.P.Val",
        stat = "t",
        label = "tf",
        nudge_factor = 0.2,
        fills = fills,
        title = title
      ),
      format = "rds"
    ),
    NULL
  ),
  tar_target(
    gsea_dot_plot,
    plot_gsea_dot(
      list(
        `TGFβ\n` = gsea_tgfb,
        `TGFβ +\nAZD` = gsea_azd,
        `TGFβ +\nVB` = gsea_vb,
        `TGFβ +\nAZD/VB` = gsea_dual
      )
    ),
    format = "rds"
  ),
  tar_target(
    edge_plot,
    plot_edge(
      list(
        TGFβ = gsea_tgfb,
        AZD = gsea_azd,
        VB = gsea_vb,
        `AZD/VB` = gsea_dual
      ),
      dds
    ),
    format = "rds"
  ),

  # mice --------------------------------------------------------------------

  tar_target(
    mice_file,
    raw_data_path("^mice.csv"),
    format = "file"
  ),
  tar_target(
    mice_raw,
    readr::read_csv(mice_file, show_col_types = FALSE)
  ),

  # mice wt -----------------------------------------------------------------

  tar_target(
    mice_wt_data,
    format_wt(mice_raw)
  ),
  tar_target(
    mice_wt_mcti_data,
    mice_wt_data |>
      dplyr::filter(
        start_date %in% c("2021-12-01", "2021-07-29", "2021-05-31", "2022-04-01")
      )
  ),
  tar_target(
    mice_wt_mcti_bin,
    bin_wt(mice_wt_mcti_data)
  ),
  tar_target(
    mice_wt_mcti_bin_stat,
    stat_bin_wt(mice_wt_mcti_bin)
  ),
  tar_target(
    mice_wt_mcti_plot,
    plot_wt_timeline(
      mice_wt_mcti_bin,
      mice_wt_mcti_bin_stat
    ),
    format = "rds"
  ),

  # mice vent ---------------------------------------------------------------

  tar_target(
    mice_vent,
    format_vent(mice_raw)
  ),
  tar_target(
    mice_vent_mcti_data,
    mice_vent |>
      dplyr::filter(
        start_date %in% c("2021-12-01", "2021-07-29", "2021-05-31", "2022-04-01")
      )
  ),
  tar_target(
    mice_vent_mcti_stats,
    stats_vent(mice_vent_mcti_data)
  ),
  tar_map(
    values = list(
      measure = c("cst", "prime3", "prime8"),
      ytitle = list(
        expression(paste("Compliance (mL/cm H" [2], "O)")),
        expression(paste("Elastance (cm H" [2], "O/mL)")),
        expression(paste("Elastance (cm H" [2], "O/mL)"))
      )
    ),
    names = measure,
    tar_target(
      mice_vent_mcti_plot,
      plot_mice(
        mice_vent_mcti_data,
        mice_vent_mcti_stats,
        measure = measure,
        ytitle = ytitle,
        y = "value_corr"
      ),
      format = "rds"
    )
  ),

  # mice histo --------------------------------------------------------------

  tar_map(
    values = list(
      col = c("Ashcroft_avg", "OHP_total"),
      names = c("ashcroft", "ohp"),
      ytitle = c("Ashcroft score", "Hydroxyproline (μg/R lung)")
    ),
    names = names,
    tar_target(
      mice,
      format_mice_col(mice_raw, col)
    ),
    tar_target(
      mice_mcti,
      dplyr::filter(mice, genotype == "WT")
    ),
    tar_target(
      mice_mcti_stats,
      stats_histo(mice_mcti, measure = names, comps = comps_bleo_2)
    ),
    tar_target(
      mice_mcti_plot,
      plot_mice(
        mice_mcti,
        mice_mcti_stats,
        ytitle = ytitle,
        measure = names,
        y = "value"
      ),
      format = "rds"
    ),
    NULL
  ),

  # mims --------------------------------------------------------------------

  tar_target(
    mims_files,
    raw_data_path("_data.csv"),
    format = "file"
  ),
  tar_target(
    mims_clean,
    clean_mims(mims_files)
  ),
  tar_map(
    values = list(
      names = c(
        "azd_1_alv_1",
        "azd_1_fib_1",
        "azd_1_fib_2",
        "azd_2_fib_1",
        "azd_3_fib_1",
        "azd_3_fib_2",
        "azd_3_fib_3",
        "azd_3_fib_4",
        "azd_3_fib_5",
        "azd_3_fib_6",
        "azd_3_fib_7",
        "azd_3_fib_8",
        "azd_3_fib_9",
        "bleo_1_alv_1",
        "bleo_1_fib_1",
        "bleo_2_fib_1",
        "bleo_2_fib_2",
        "bleo_3_fib_1",
        "bleo_3_fib_2",
        "bleo_3_fib_3",
        "bleo_3_fib_4",
        "bleo_3_fib_5",
        "bleo_3_fib_6",
        "bleo_3_fib_7",
        "control_1_alv_1",
        "control_2_alv_1",
        "control_2_alv_2",
        "control_3_alv_1",
        "control_3_alv_2",
        "control_3_alv_3",
        "vb_1_alv_1",
        "vb_1_fib_1",
        "vb_2_fib_1",
        "vb_3_fib_1",
        "vb_3_fib_2",
        "vb_3_fib_3",
        "vb_3_fib_4",
        "vb_3_fib_5",
        "vb_3_fib_6",
        "vb_3_fib_7",
        "vb_3_fib_8",
        "vb_3_fib_9"
      )
    ),
    names = names,
    tar_target(
      mims_files,
      raw_data_path(stringr::str_c(names, ".+(\\d|r)\\.nrrd"))
    ),
    tar_target(
      mims_img,
      read_mims(mims_files)
    ),
    tar_target(
      mims_ratios,
      extract_ratios(mims_img, names) |>
        format_ratios()
    ),
    tar_target(
      mims_averages,
      average_mims(mims_ratios)
    ),
    NULL
  ),
  tar_target(
    mims_averages,
    dplyr::bind_rows(
      mims_averages_azd_1_alv_1,
      mims_averages_azd_1_fib_1,
      mims_averages_azd_1_fib_2,
      mims_averages_azd_2_fib_1,
      mims_averages_azd_3_fib_1,
      mims_averages_azd_3_fib_2,
      mims_averages_azd_3_fib_3,
      mims_averages_azd_3_fib_4,
      mims_averages_azd_3_fib_5,
      mims_averages_azd_3_fib_6,
      mims_averages_azd_3_fib_7,
      mims_averages_azd_3_fib_8,
      mims_averages_azd_3_fib_9,
      mims_averages_bleo_1_alv_1,
      mims_averages_bleo_1_fib_1,
      mims_averages_bleo_2_fib_1,
      mims_averages_bleo_2_fib_2,
      mims_averages_bleo_3_fib_1,
      mims_averages_bleo_3_fib_2,
      mims_averages_bleo_3_fib_3,
      mims_averages_bleo_3_fib_4,
      mims_averages_bleo_3_fib_5,
      mims_averages_bleo_3_fib_6,
      mims_averages_bleo_3_fib_7,
      mims_averages_control_1_alv_1,
      mims_averages_control_2_alv_1,
      mims_averages_control_2_alv_2,
      mims_averages_control_3_alv_1,
      mims_averages_control_3_alv_2,
      mims_averages_control_3_alv_3,
      mims_averages_vb_1_alv_1,
      mims_averages_vb_1_fib_1,
      mims_averages_vb_2_fib_1,
      mims_averages_vb_3_fib_1,
      mims_averages_vb_3_fib_2,
      mims_averages_vb_3_fib_3,
      mims_averages_vb_3_fib_4,
      mims_averages_vb_3_fib_5,
      mims_averages_vb_3_fib_6,
      mims_averages_vb_3_fib_7,
      mims_averages_vb_3_fib_8,
      mims_averages_vb_3_fib_9
    )
  ),
  tar_target(
    mims_summary,
    summarize_mims(mims_averages)
  ),
  tar_target(
    mims_stats,
    stats_mims(mims_summary)
  ),
  tar_target(
    mims_plot,
    plot_mice(
      mims_summary,
      stats = mims_stats,
      measure = "tracing",
      x = "group",
      y = "value_corr",
      ytitle = "Isotope enrichment\n(normalized)"
    ) +
      ggplot2::facet_wrap(
        facets = ggplot2::vars(tracer),
        scales = "free_y",
        nrow = 1
      ),
    format = "rds"
  ),

  # VB253 -------------------------------------------------------------------

  tar_target(
    vb_invitro_sma_file,
    raw_data_path("invitro_sma.csv"),
    format = "file"
  ),
  tar_target(
    vb_invitro_sma_data,
    read_vb_invitro_sma(vb_invitro_sma_file)
  ),
  tar_map(
    values = list(
      type = c("IPF", "Ctl"),
      names = c("ipf", "ctl")
    ),
    names = names,
    tar_map(
      values = list(
        y = c("sma_norm", "nuclei_norm"),
        ytitle = c("α-SMA fold change", "Nuclei fold change"),
        names = c("sma", "nuclei")
      ),
      names = names,
      tar_target(
        vb_invitro_plot,
        plot_vb_invitro_sma(vb_invitro_sma_data, type, y, ytitle),
        format = "rds"
      )
    )
  ),
  tar_target(
    vb_invitro_smad_file,
    raw_data_path("invitro_smad3.csv"),
    format = "file"
  ),
  tar_target(
    vb_invitro_smad_data,
    read_vb_invitro_sma(vb_invitro_smad_file)
  ),
  tar_target(
    vb_invitro_smad_plot,
    plot_vb_invitro_sma(vb_invitro_smad_data, "IPF", ytitle = "Nuclear Smad3"),
    format = "rds"
  ),
  tar_target(
    vb_mice_file,
    raw_data_path("vb_mice.csv"),
    format = "file"
  ),
  tar_target(
    vb_mice_data,
    read_vb_mice(vb_mice_file)
  ),
  tar_map(
    values =
      tibble::tribble(
        ~ages,   ~measure,   ~ytitle,                       ~comps,
        "young", "ashcroft", "Ashcroft score",              rlang::sym("comps_bleo_3"),
        "young", "sma",      "α-SMA area",                  rlang::sym("comps_bleo_3"),
        "young", "Penh",     "Penh",                        rlang::sym("comps_bleo_4"),
        "old",   "ashcroft", "Ashcroft score",              rlang::sym("comps_bleo_4"),
        "old",   "sma",      "α-SMA area",                  rlang::sym("comps_bleo_4"),
        "old",   "lactate",  "Lactate (nmol / mg protein)", rlang::sym("comps_bleo_4")
      ) |>
      tidyr::unite(names, ages, measure, remove = FALSE),
    names = names,
    tar_target(
      vb_mice,
      dplyr::filter(vb_mice_data, age == ages, measurement == measure)
    ),
    tar_target(
      vb_mice_stats,
      stats_histo(vb_mice, measure = names, comps = comps) |>
        dplyr::mutate(measurement = measure)
    ),
    tar_target(
      vb_mice_plot,
      plot_mice(
        vb_mice,
        vb_mice_stats,
        ytitle = ytitle,
        measure = measure,
        y = "value"
      ),
      format = "rds"
    ),
    NULL
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
  tar_quarto(
    redox_analysis,
    path = "analysis/redox.qmd"
  ),
  tar_quarto(
    bleo_mid_analysis,
    path = "analysis/bleo-mid.qmd"
  ),
  tar_quarto(
    mims_analysis,
    path = "analysis/mims.qmd"
  ),
  tar_quarto(
    bcecf_analysis,
    path = "analysis/bcecf.qmd"
  ),

  # figure images -----------------------------------------------------------

  tar_map(
    values = tibble::tribble(
      ~path,                        ~scale, ~hjust, ~vjust, ~names,
      "ipf-mct-blots.png",          1.1,    0,      -0.05,  "mct_ipf",
      "bleo-mct-blots.png",         1.1,    0,      -0.05,  "mct_bleo",
      "tgfb-mct-blots.png",         1.1,    0,      -0.05,  "mct_tgfb",
      "sirna-mct-blots.png",        1,      0,      0,      "sma_sirna",
      "dual-azd-ipf-lf-blot-2.png", 0.9,    0,      0,      "ipf_lf_azd",
      "dual-azd-sma-blot.png",      1,      0,      -0.05,  "sma_azd",
      "dual-ar-sma-blot.png",       1.2,    0.1,    -0.05,  "sma_ar",
      "dual-azd-contract.png",      1,      0,      -0.05,  "contraction",
      "bleo-trichrome.png",         1.1,    0,      0,      "trichrome",
      "mct-bleo-timeline.png",          1,      0,      0,  "timeline",
      "dual-azd-smad3-blot.png",    1,      0,      -0.05,  "smad3_azd",
      "dual-azd-erk-blot.png",      1,      0,      -0.05,  "erk_azd",
      "dual-azd-hif1a-blot.png",    1,      0,      -0.05,  "hif_azd",
      "lactate-sma-blot.png",       1.2,    0,      -0.05,  "sma_lac",
      "dual-azd-kla-h3.png",        1,      0,      0,      "azd_kla",
      "mims-panel.png",             1.05,   0.05,   0,      "mims",
      "vb-bleo-timeline.png",       1,      0,      0,      "vb_timeline",
      "bleo-vb-young-sma.png",      1,      0,      0,      "vb_young_sma",
      "bleo-vb-aged-sma.png",       1,      0,      0,      "vb_aged_sma",
      "glucose-tracer.png",         1,      0,      0,      "tracer"
    ),
    names = names,
    tar_target(
      fig_img_file,
      manuscript_path(path),
      format = "file"
    ),
    tar_target(
      fig_img,
      plot_image(fig_img_file),
      format = "rds"
    ),
    NULL
  ),

  # figure 1 ----------------------------------------------------------------

  tar_target(
    fig01,
    make_fig01(
      fig_img_mct_ipf,
      blot_plot_fig_ipf,
      fig_img_mct_bleo,
      blot_plot_fig_bleo,
      fig_img_mct_tgfb,
      blot_plot_fig_tgfb
    ) |>
      write_figures("F_01"),
    format = "file"
  ),

  # figure 2 ----------------------------------------------------------------

  tar_target(
    fig02,
    make_fig02(
      fig_img_sma_sirna,
      blot_plot_fig_sirna,
      fig_img_ipf_lf_azd,
      blot_plot_fig_ipf_lf,
      fig_img_sma_azd,
      blot_plot_fig_azd +
        ggplot2::guides(fill = "none"),
      fig_img_contraction,
      contraction_plot,
      rnaseq_pca_plot,
      gsea_dot_plot,
      deg_vol_dual,
      edge_plot
    ) |>
      write_figures("F_02"),
    format = "file"
  ),

  # figure 3 ----------------------------------------------------------------

  tar_target(
    fig03,
    make_fig03(
      plate_plot_lactate_sirna,
      plate_plot_lactate_azd,
      seahorse_mcti_timeline_plot,
      seahorse_mcti_atp_bars,
      seahorse_mcti_atp_pheno
    ) |>
      write_figures("F_03"),
    format = "file"
  ),

  # figure 4 ----------------------------------------------------------------

  tar_target(
    fig04,
    make_fig04(
      metab_mcti_tar_pca_filter_extra +
        ggplot2::theme(legend.margin = ggplot2::margin(t = 0)),
      metab_mcti_tar_lactate_extra,
      metab_mcti_tar_vol_dual_extra,
      msea_table_img_dual_extra,
      metab_mcti_tar_pca_filter_intra +
        ggplot2::theme(legend.margin = ggplot2::margin(t = -40)),
      metab_mcti_tar_lactate_intra,
      metab_mcti_tar_vol_dual_intra,
      msea_table_img_dual_intra,
      fig_img_tracer,
      mid_mcti_moi_plot_glc6,
      mid_mcti_moi_plot_lac3,
      nad_plot_nad_norm,
      nad_plot_nadh_norm,
      nad_plot_nadh_ratio,
      ros_plot_cellrox,
      ros_plot_ratio
    ) |>
      write_figures("F_04"),
    format = "file"
  ),

  # figure 5 ----------------------------------------------------------------

  tar_target(
    fig05,
    make_fig05(
      fig_img_smad3_azd,
      blot_plot_fig_psmad3,
      fig_img_erk_azd,
      blot_plot_fig_perk
    ) |>
      write_figures("F_05"),
    format = "file"
  ),

  # figure 6 ----------------------------------------------------------------

  tar_target(
    fig06,
    make_fig06(
      fig_img_timeline,
      mice_vent_mcti_plot_cst,
      mice_vent_mcti_plot_prime8,
      fig_img_trichrome,
      mice_mcti_plot_ashcroft,
      mice_mcti_plot_ohp
    ) |>
      write_figures("F_06"),
    format = "file"
  ),

  # figure 7 ----------------------------------------------------------------

  tar_target(
    fig07,
    make_fig07(
      metab_bleo_lacate_plot,
      fig_img_mims,
      mims_plot
    ) |>
      write_figures("F_07"),
    format = "file"
  ),

  # figure 8 ----------------------------------------------------------------

  tar_target(
    fig08,
    make_fig08(
      vb_invitro_plot_sma_ipf,
      vb_invitro_plot_nuclei_ipf,
      fig_img_vb_timeline,
      vb_mice_plot_young_Penh,
      vb_mice_plot_young_ashcroft,
      fig_img_vb_young_sma,
      vb_mice_plot_young_sma +
        ggplot2::scale_y_continuous(
          labels = scales::label_percent(),
          breaks = scales::breaks_extended(n = 6, only.loose = TRUE),
          expand = ggplot2::expansion(c(0, 0)),
          limits = nice_limits
        ),
      vb_mice_plot_old_ashcroft,
      fig_img_vb_aged_sma,
      vb_mice_plot_old_sma +
        ggplot2::scale_y_continuous(
          labels = scales::label_percent(),
          breaks = scales::breaks_extended(n = 6, only.loose = TRUE),
          expand = ggplot2::expansion(c(0, 0)),
          limits = nice_limits
        ),
      vb_mice_plot_old_lactate
    ) |>
      write_figures("F_08"),
    format = "file"
  ),

  # figure 1S ---------------------------------------------------------------

  tar_target(
    fig01s,
    make_fig01s(
      pg_plot_sirna,
      fig_img_sma_ar,
      blot_plot_fig_ar +
        ggplot2::guides(fill = "none"),
      pg_plot_azd,
      pg_plot_ar
    ) |>
      write_figures("SF_01"),
    format = "file"
  ),

  # figure 2S ---------------------------------------------------------------

  tar_target(
    fig02s,
    make_fig02s(
      deg_vol_tgfb,
      deg_vol_azd,
      deg_vol_vb
    ) |>
      write_figures("SF_02"),
    format = "file"
  ),

  # figure 3S ---------------------------------------------------------------

  tar_target(
    fig03s,
    make_fig03s(
      plate_plot_lactate_ar,
      plate_plot_glucose_azd,
      seahorse_mcti_summary_plot,
      seahorse_mcti_stress_plot
    ) |>
      write_figures("SF_03"),
    format = "file"
  ),

  # figure 4S ---------------------------------------------------------------

  tar_target(
    fig04s,
    make_fig04s(
      metab_mcti_tar_vol_tgfb_extra,
      msea_table_img_tgfb_extra,
      metab_mcti_tar_vol_tgfb_intra,
      msea_table_img_tgfb_intra,
      metab_mcti_tar_vol_azd_extra,
      msea_table_img_azd_extra,
      metab_mcti_tar_vol_azd_intra,
      msea_table_img_azd_intra,
      metab_mcti_tar_vol_vb_extra,
      msea_table_img_vb_extra,
      metab_mcti_tar_vol_vb_intra,
      msea_table_img_vb_intra,
      mid_mcti_moi_plot_gln5
    ) |>
      write_figures("SF_04"),
    format = "file"
  ),

  # figure 6S ---------------------------------------------------------------

  tar_target(
    fig06s,
    make_fig06s(
      nad_plot_nadp_norm,
      nad_plot_nadph_norm,
      nad_plot_nadph_ratio,
      ros_plot_mitotracker,
      metab_mcti_tar_proline_intra,
      mid_plot_gln5_proline
    ) |>
      write_figures("SF_06"),
    format = "file"
  ),

  # figure 7S ---------------------------------------------------------------

  tar_target(
    fig07s,
    make_fig07s(
      fig_img_sma_lac,
      blot_plot_lactate_sma,
      fig_img_hif_azd,
      blot_plot_six_hif1a
    ) |>
      write_figures("SF_07"),
    format = "file"
  ),

  # figure 8S ---------------------------------------------------------------

  tar_target(
    fig08s,
    make_fig08s(
      mice_wt_mcti_plot
    ) |>
      write_figures("SF_08"),
    format = "file"
  ),

  # figure 9S ---------------------------------------------------------------

  tar_target(
    fig09s,
    make_fig09s(
      metab_bleo_tar_vol_bleo_plasma,
      msea_table_img_bleo_plasma,
      metab_bleo_tar_vol_bleo_lung,
      msea_table_img_bleo_lung,
      metab_bleo_tar_vol_azd_plasma,
      msea_table_img_azd_plasma,
      metab_bleo_tar_vol_azd_lung,
      msea_table_img_azd_lung,
      metab_bleo_tar_vol_vb_plasma,
      msea_table_img_vb_plasma,
      metab_bleo_tar_vol_vb_lung,
      msea_table_img_vb_lung
    ) |>
      write_figures("SF_09"),
    format = "file"
  ),
  # tar_target(
  #   fig09s2,
  #   make_fig09s2(
  #     mid_bleo_plasma_glucose_plot,
  #     mid_bleo_plasma_ratio
  #   ) |>
  #     write_figures("Figure 09.supplement.2"),
  #   format = "file"
  # ),

  # figure 10S --------------------------------------------------------------

  tar_target(
    fig10s,
    make_fig10s(
      vb_invitro_smad_plot
    ) |>
      write_figures("SF_10"),
    format = "file"
  ),

  # resources ---------------------------------------------------------------

  tar_target(
    resources_file,
    manuscript_path("resources.csv"),
    format = "file"
  ),
  tar_target(
    resources_table,
    create_resources(resources_file)
  ),

  # manuscript --------------------------------------------------------------

  tar_target(
    template,
    manuscript_path("template.docx"),
    format = "file"
  ),
  tar_target(
    pkg_citations,
    write_pkg_cites(),
    cue = tar_cue(mode = "always")
  ),
  tar_target(
    csl,
    manuscript_path("jci.csl"),
    format = "file"
  ),
  tar_target(
    refs_ms,
    rbbt::bbt_update_bib(
      path = "manuscript/manuscript.qmd",
      ignore = c("R-base"),
      path_bib = "manuscript/ms.bib"
    ),
    cue = tar_cue("always")
  ),
  tar_target(
    refs_supp,
    rbbt::bbt_update_bib(
      path = "manuscript/supplement.qmd",
      ignore = c("R-base"),
      path_bib = "manuscript/supp.bib"
    ),
    cue = tar_cue("always")
  ),
  tar_quarto(
    manuscript,
    path = manuscript_path("manuscript.qmd"),
  ),
  tar_quarto(
    supplement,
    path = manuscript_path("supplement.qmd"),
  ),

  NULL
)
