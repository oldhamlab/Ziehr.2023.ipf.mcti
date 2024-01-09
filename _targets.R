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
    "scales",
    "SummarizedExperiment"
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

  # nad ---------------------------------------------------------------------

  tar_target(
    nad_clean,
    clean_nad(plates_conc_nad_protein, plates_conc_nad_nad)
  ),
  tar_map(
    values = list(
      names = c(
        "nad",
        "nad_norm",
        "nadh",
        "nadh_norm",
        "ratio"
      ),
      measure = c(
        "NAD",
        "NAD_norm",
        "NADH",
        "NADH_norm",
        "Ratio"
      ),
      ytitle = c(
        "NAD (pmol)",
        "NAD (pmol / μg protein)",
        "NADH (pmol)",
        "NADH (pmol / μg protein)",
        "NADH / NAD"
      )
    ),
    names = names,
    tar_target(
      nad_filter,
      dplyr::filter(nad_clean, measurement == measure)
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
      ),
      format = "rds"
    )
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

  tar_target(
    metab_mcti_tar_file,
    raw_data_path("myofib_azd3965-vb124_2022-02-10.xlsx"),
    format = "file"
  ),
  tar_target(
    metab_mcti_tar_raw,
    format_metab_targeted(metab_mcti_tar_file)
  ),
  tar_target(
    metab_mcti_samples,
    make_mcti_samples()
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
    metab_mcti_tar_clean[, metab_mcti_tar_clean$replicate != 5]
  ),
  tar_target(
    metab_mcti_tar_treated,
    metab_mcti_tar_out[, metab_mcti_tar_out$treatment != "None"]
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
        c("Ctl", "TGFβ"),
        c("Veh", "AZD"),
        c("Veh", "VB"),
        c("Veh", "AZD/VB")
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
      plot_vol_metab(
        metab_mcti_tar_tt,
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
      metab_mcti_tar_msea_plot,
      plot_msea(metab_mcti_tar_msea, title = title, colors = fills),
      format = "rds"
    ),
    tar_target(
      metab_mcti_tar_msea_table,
      plot_msea_table(metab_mcti_tar_msea, title = title, clr = fills),
      format = "rds"
    ),
    tar_target(
      msea_table_file,
      write_table(
        metab_mcti_tar_msea_table,
        path = "analysis/figures/msea/",
        filename
      ),
      format = "file"
    ),
    tar_target(
      msea_table_img,
      plot_table(msea_table_file),
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
  tar_target(
    metab_mcti_tar_lactate,
    plot_mois(metab_mcti_tar_treated, metab_mcti_tar_tt, "lactate"),
    format = "rds"
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
      names = c("lung", "plasma"),
      input = rlang::syms(c(
        "metab_bleo_lung_tar_file",
        "metab_bleo_plasma_tar_file"
      ))
    ),
    names = names,
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
          c("Ctl", "Bleo"),
          c("Bleo", "AZD"),
          c("Bleo", "VB")
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
        plot_vol_metab(metab_bleo_tar_tt, fills, title),
        format = "rds"
      ),
      tar_target(
        metab_bleo_tar_msea,
        run_msea(metab_bleo_tar_tt, metabolite_pathways)
      ),
      tar_target(
        bleo_tar_msea_plot,
        plot_msea(metab_bleo_tar_msea, colors = unname(fills), src = "KEGG"),
        format = "rds"
      ),
      tar_target(
        metab_bleo_tar_msea_table,
        plot_msea_table(metab_bleo_tar_msea, title = title, clr = fills),
        format = "rds"
      ),
      tar_target(
        msea_table_file,
        write_table(
          metab_bleo_tar_msea_table,
          path = "analysis/figures/msea/",
          filename
        ),
        format = "file"
      ),
      tar_target(
        msea_table_img,
        plot_table(msea_table_file),
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

  # rnaseq ------------------------------------------------------------------

  tar_target(
    dds,
    count_rnaseq()
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
    decoupler_matrix,
    prep_counts(dds)
  ),
  tar_target(
    progeny_network,
    decoupleR::get_progeny(organism = "human")
  ),
  tar_target(
    progeny_results,
    decoupleR::run_mlm(
      mat = decoupler_matrix,
      net = progeny_network,
      .mor = "weight"
    )
  ),
  tar_target(
    progeny_plot,
    plot_progeny(progeny_results),
    format = "rds"
  ),
  tar_target(
    tf_network,
    decoupleR::get_collectri(organism = "human")
  ),
  tar_target(
    tf_results,
    run_tfea(decoupler_matrix, tf_network)
  ),
  tar_target(
    tf_plot,
    plot_tf(tf_results),
    format = "rds"
  ),
  tar_map(
    values = list(
      names = c(
        "main_tgfb",
        "main_azd",
        "main_vb",
        "main_dual",
        "int_azd",
        "int_vb",
        "int_dual",
        "tgfb_azd",
        "tgfb_vb",
        "tgfb_dual"
      ),
      con = list(
        c("condition", "TGFβ", "Ctl"),
        c("treatment", "AZD", "Veh"),
        c("treatment", "VB", "Veh"),
        c("treatment", "AZD/VB", "Veh"),
        list("conditionTGFβ.treatmentAZD"),
        list("conditionTGFβ.treatmentVB"),
        list("conditionTGFβ.treatmentAZD.VB"),
        list(c("treatment_AZD_vs_Veh", "conditionTGFβ.treatmentAZD")),
        list(c("treatment_VB_vs_Veh", "conditionTGFβ.treatmentVB")),
        list(c("treatment_AZD.VB_vs_Veh", "conditionTGFβ.treatmentAZD.VB"))
      ),
      title = c(
        "TGFβ v. Ctl",
        "AZD v. Veh",
        "VB v. Veh",
        "AZD/VB v. Veh",
        "AZD × TGFβ",
        "VB × TGFβ",
        "AZD/VB × TGFβ",
        "AZD v. Veh in TGFβ",
        "VB v. Veh in TGFβ",
        "AZD/VB v. Veh in TGFβ"
      ),
      bins = list(
        c(30, 10),
        c(20, 80),
        c(45, 15),
        c(45, 15),
        c(75, 2000),
        c(100, 60),
        c(100, 60),
        c(30, 110),
        c(30, 15),
        c(30, 10)
      ),
      fills = list(
        c("Ctl", "TGFβ"),
        c("Veh", "AZD"),
        c("Veh", "VB"),
        c("Veh", "AZD/VB"),
        c("Veh", "AZD"),
        c("Veh", "VB"),
        c("Veh", "AZD/VB"),
        c("Veh", "AZD"),
        c("Veh", "VB"),
        c("Veh", "AZD/VB")
      )
    ),
    names = names,
    tar_target(
      deg,
      rnaseq_results(dds, con = con)
    ),
    tar_target(
      deg_vol,
      plot_vol_rnaseq(deg, title = title, bins = bins),
      format = "rds"
    ),
    tar_target(
      gsea,
      run_gsea(deg, gsea_pathways)
    ),
    tar_target(
      gsea_table,
      plot_gsea_table(gsea, "HALLMARK", title = title, clr = fills),
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
      plot_table(gsea_table_file),
      format = "rds"
    ),
    tar_target(
      tfea,
      score_tf(deg, tf_network)
    ),
    tar_target(
      tfea_vol,
      plot_vol_tfea(tfea, fills, title),
      format = "rds"
    ),
    NULL
  ),
  tar_target(
    gsea_dot_plot,
    plot_gsea_dot(
      list(
        `TGFβ\n` = gsea_main_tgfb,
        `TGFβ +\nAZD` = gsea_tgfb_azd,
        `TGFβ +\nVB` = gsea_tgfb_vb,
        `TGFβ +\nAZD/VB` = gsea_tgfb_dual
      )
    ),
    format = "rds"
  ),
  tar_target(
    edge_plot,
    plot_edge(
      list(
        TGFβ = gsea_main_tgfb,
        AZD = gsea_tgfb_azd,
        VB = gsea_tgfb_vb,
        `AZD/VB` = gsea_tgfb_dual
      ),
      dds
    ),
    format = "rds"
  ),

  # mice --------------------------------------------------------------------

  tar_target(
    mice_file,
    raw_data_path("mice.csv"),
    format = "file"
  ),
  tar_target(
    mice_raw,
    readr::read_csv(mice_file, show_col_types = FALSE)
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

  NULL
)
