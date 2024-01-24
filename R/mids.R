# mid.R

# functions ---------------------------------------------------------------

isotope_library <-
  tibble::tribble(
    ~ metabolite, ~ formula,   ~ polarity,
    "2HG",        "C5H8O5",     "negative",
    "2OG",        "C5H6O5",     "negative",
    "ALA",        "C3H7NO2",    "positive",
    "ASP",        "C4H7NO4",    "negative",
    "CIT",        "C6H8O7",     "negative",
    "GLU",        "C5H9NO4",    "positive",
    "GLN",        "C5H10N2O3",  "positive",
    "LAC",        "C3H6O3",     "negative",
    "MAL",        "C4H6O5",     "negative",
    "PYR",        "C3H4O3",     "negative",
    "SER",        "C3H7NO3",    "positive",
    "SUC",        "C4H6O4",     "negative",
    "3PG",        "C3H7O7P",    "negative",
    "ACO",        "C6H6O6",     "negative",
    "FBP",        "C6H14O12P2", "negative",
    "6PG",        "C6H13O10P",  "negative",
    "GLY",        "C2H5NO2",    "positive",
    "PRO",        "C5H9NO2",    "positive",
    "G3P",        "C3H9O6P",    "negative",
    "GLC",        "C6H12O6",    "negative"
    # "palmitate", "C16H32O2", "negative",
    # "PEP", "C3H5O6P", "negative",
    # "sedoheptulose", "C7H14O7", "negative",
    # "DHAP", "C3H7O6P", "negative",
    # "GAP", "C3H7O6P", "negative",
    # "G1P", "C6H13O9P", "negative",
    # "G6P", "C6H13O9P", "negative",
    # "R5P", "C5H11O8P", "negative"
  )

make_correction_matrices <- function(isotope_library) {
  isotope_library |>
    dplyr::mutate(
      matrix = purrr::map2(formula, polarity, mzrtools::mz_iso_quant),
      matrix = purrr::map(matrix, purrr::pluck, "prob_matrix")
    )
}

mmult <- function(m, df) {
  mid <- df$mid
  if (nrow(m) > length(mid)) {
    m <- m[1:length(mid), 1:length(mid)]
  }
  mid_corr <- mzrtools::mz_iso_correct(m, mid)
  dplyr::bind_cols(df, mid_corr = mid_corr)
}

correct_mids <- function(mid, correction_matrices) {
  mid |>
    tidyr::nest() |>
    dplyr::inner_join(
      dplyr::select(correction_matrices, -c("formula", "polarity")),
      by = "metabolite"
    ) |>
    dplyr::mutate(data = purrr::map2(.data$matrix, .data$data, mmult)) |>
    tidyr::unnest(c("data")) |>
    dplyr::filter(.data$isotope %in% stringr::str_c("M", 0:6)) |>
    dplyr::select(-"matrix")
}

plot_mids_all <- function(df, label) {
  title <-
    list(
      glc2 = expression(paste("[1,2-"^13, "C"[2], "]-glucose")),
      gln5 = expression(paste("[U-"^13, "C"[5], "]-glutamine")),
      glc6 = expression(paste("[U-"^13, "C"[6], "]-glucose")),
      lac3 = expression(paste("[U-"^13, "C"[3], "]-lactate"))
    )

  x <-
    df |>
    dplyr::filter(tracer %in% label) |>
    dplyr::group_by(condition, treatment, metabolite, isotope) |>
    dplyr::summarize(
      mean = mean(mid_corr),
      se = sd(mid) / sqrt(dplyr::n())
    ) |>
    dplyr::group_by(condition, treatment, metabolite) |>
    dplyr::mutate(
      means = 1 - cumsum(mean) + mean,
      ymin = means + se,
      ymax = means - se,
      group = stringr::str_c(condition, treatment, sep = "."),
      group = factor(group, levels = c(
        "Ctl.Veh", "Ctl.AZD", "Ctl.VB", "Ctl.AZD/VB",
        "TGFβ.Veh", "TGFβ.AZD", "TGFβ.VB", "TGFβ.AZD/VB"
      ))
    )

  ggplot2::ggplot(x) +
    ggplot2::facet_wrap(~ metabolite) +
    ggplot2::aes(
      x = group,
      y = mean,
      fill = isotope
    ) +
    ggplot2::geom_col() +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = ymin,
        ymax = ymax
      ),
      position = ggplot2::position_dodge(width = 0.75),
      color = "grey50",
      width = 0,
      linewidth = 0.25
    ) +
    ggplot2::scale_fill_viridis_d() +
    ggplot2::scale_y_continuous(
      breaks = seq(0, 1, 0.2),
      expand = ggplot2::expansion(c(0, 0))
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Mole fraction",
      fill = NULL,
      title = title[[label]]
    ) +
    theme_plot() +
    ggplot2::theme(
      legend.position = "right",
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)
    )
}

# qbias -------------------------------------------------------------------

format_qbias <- function(file_name){
  readr::read_csv(
    file_name,
    show_col_types = FALSE
  ) |>
    dplyr::rename(ion = ...1) |>
    tidyr::pivot_longer(!ion, names_to = "filename", values_to = "area") |>
    dplyr::filter(!is.na(area)) |>
    dplyr::mutate(
      window = as.integer(stringr::str_extract(filename, "(?<=W)\\d{1}")),
      replicate = as.integer(stringr::str_extract(filename, "(?<=_)\\d{1}"))
    ) |>
    dplyr::select(-filename) |>
    tidyr::separate(
      .data$ion,
      into = c("metabolite", "isotope"),
      sep = " ",
      convert = TRUE
    ) |>
    dplyr::mutate(carbons = dplyr::case_when(
      .data$metabolite %in% c("CIT", "FBP", "6PG") ~ 6,
      .data$metabolite %in% c("2HG", "2OG", "GLU", "GLN", "PRO") ~ 5,
      .data$metabolite %in% c("ASP", "MAL", "SUC") ~ 4,
      .data$metabolite %in% c("LAC", "PYR", "ALA", "SER") ~ 3,
      .data$metabolite %in% c("GLY") ~ 2
    )) |>
    dplyr::filter(.data$window <= .data$carbons) |>
    tidyr::pivot_wider(names_from = .data$isotope, values_from = .data$area) |>
    dplyr::mutate(ratio = .data$M1 / .data$M0) |>
    dplyr::group_by(.data$metabolite)
}

calculate_predicted_ratios <- function(isotope_library){
  isotope_library |>
    dplyr::mutate(table = purrr::map2(
      formula,
      polarity,
      \(x, y) mzrtools::mz_iso_quant(molecule = x, polarity = y)[["prob_matrix"]]
    )) |>
    dplyr::mutate(pred_ratio = purrr::map_dbl(table, \(x) x[[2, 1]]/x[[1, 1]])) |>
    dplyr::select(.data$metabolite, .data$pred_ratio)
}

calculate_correction_factors <- function(qbias_ratios, pred_ratios){
  qbias_ratios |>
    dplyr::group_by(metabolite) |>
    tidyr::nest() |>
    dplyr::mutate(
      model = purrr::map(
        .data$data,
        \(x) MASS::rlm(ratio ~ poly(window, 2), data = x)
      ),
      predict = purrr::map2(.data$model, .data$data, predict)
    ) |>
    tidyr::unnest(c(data, predict)) |>
    dplyr::select("metabolite", "window", "carbons", "predict") |>
    dplyr::distinct() |>
    dplyr::left_join(pred_ratios, by = "metabolite") |>
    dplyr::mutate(cf = .data$predict / .data$pred_ratio) |>
    dplyr::select("metabolite", M = "window", "carbons", "cf") |>
    dplyr::filter(.data$M < .data$carbons) |>
    dplyr::ungroup() |>
    dplyr::mutate(M = M + 1) |>
    tidyr::pivot_wider(names_from = .data$M, values_from = .data$cf) |>
    dplyr:: mutate(
      M0 = 1,
      M1 = 1 / .data$`1` * .data$M0,
      M2 = 1 / .data$`2` * .data$M1,
      M3 = 1 / .data$`3` * .data$M2,
      M4 = 1 / .data$`4` * .data$M3,
      M5 = 1 / .data$`5` * .data$M4,
      M6 = 1 / .data$`6` * .data$M5
    ) |>
    dplyr::select("metabolite", tidyselect::matches("M[0-9]+")) |>
    tidyr::pivot_longer(
      cols = matches("M[0-9]+"),
      names_to = "isotope",
      values_to = "corr_fac",
      values_drop_na = TRUE
    ) |>
    dplyr::arrange(.data$metabolite)
}

make_correction_matrices <- function(isotope_library) {
  isotope_library |>
    dplyr::mutate(
      matrix = purrr::map2(formula, polarity, mzrtools::mz_iso_quant),
      matrix = purrr::map(matrix, purrr::pluck, "prob_matrix")
    )
}

# analysis ----------------------------------------------------------------

format_mid_mcti <- function(files, correction_factors){
  nms <- stringr::str_replace(basename(files), "\\.xlsx", "")
  files <- rlang::set_names(files, nms)

  sample_no <- 1:32
  replicate <- rep(1:4, each = 8)
  condition <- rep(c("tgfb", "control"), each = 4)
  treatment <- rep(c("DMSO", "VB124", "AZD3965", "AZD/VB"), 2)
  samples <-
    tibble::tibble(
      sample = sample_no,
      replicate = replicate,
      condition = rep(condition, 4),
      treatment = rep(treatment, 4)
    ) |>
    refactor()

  purrr::map_dfr(
    files,
    readxl::read_excel,
    sheet = 2,
    .id = "file"
  ) |>
    dplyr::select(
      file,
      id = `Raw File Name`,
      metabolite = `Compound Name`,
      area = `Peak Area`
    ) |>
    tidyr::separate(id, c("type", "sample"), sep = "_", convert = TRUE) |>
    tidyr::separate(metabolite, c("metabolite", "isotope"), sep = " ") |>
    dplyr::filter(type != "blank") |>
    tidyr::separate(file, c("run_date", "remove", "tracer"), sep = "_") |>
    dplyr::select(-c(type, remove)) |>
    dplyr::filter(!(metabolite %in% c("FBP", "6PG") & tracer == "gln5")) |>
    dplyr::left_join(samples, by = "sample") |>
    dplyr::relocate(replicate, condition, treatment, .before = "metabolite") |>
    dplyr::group_by(run_date, tracer, replicate, sample, metabolite) |>
    dplyr::left_join(correction_factors, by = c("metabolite", "isotope")) |>
    dplyr::mutate(
      area_corr = area * corr_fac,
      mid = area_corr / sum(area_corr)
    ) |>
    dplyr::group_by(dplyr::across("run_date":"metabolite"))
}

select_mid_mois <- function(df, trace) {
  metabs <- c("PYR", "LAC", "CIT", "2OG", "SUC", "MAL")

  df |>
    dplyr::filter(.data$tracer == trace) |>
    dplyr::filter(.data$metabolite %in% metabs) |>
    dplyr::filter(.data$isotope == "M0") |>
    dplyr::group_by(.data$metabolite) |>
    dplyr::mutate(grand_mean = mean(.data$mid_corr)) |>
    dplyr::group_by(.data$metabolite, .data$replicate) |>
    dplyr::mutate(
      exp_mean = mean(.data$mid_corr),
      cf = grand_mean - exp_mean,
      # labeled = 1 - .data$mid_corr + .data$cf,
      labeled = 1 - .data$mid_corr,
      metabolite = factor(metabolite, levels = metabs)
    ) |>
    dplyr::select(
      "sample",
      group = "replicate",
      "condition":"metabolite",
      "area",
      "mid_corr",
      "labeled"
    ) |>
    dplyr::group_by(.data$metabolite)
}


format_bleo_mids <- function(files, mice) {
  nms <- stringr::str_extract(basename(files), "lung|plasma")
  files <- rlang::set_names(files, nms)

  sample <- as.character(1:17)
  start_date <- lubridate::as_date(rep("2022-08-05", 17))
  mouse_no <- c(1:6, 8, 10:11, 13:18, 23, 24)

  samples <-
    tibble::tibble(
      sample = sample,
      start_date = start_date,
      mouse_no = mouse_no
    ) |>
    dplyr::left_join(mice, by = c("start_date", "mouse_no")) |>
    dplyr::mutate(
      group = factor(
        group,
        levels = c("Veh", "Bleo", "AZD", "VB"),
        labels = c("Ctl", "Bleo", "AZD", "VB")
      )
    )

  purrr::map_dfr(
    files,
    readxl::read_excel,
    sheet = 2,
    .id = "tissue"
  ) |>
    dplyr::select(
      tissue,
      id = `Raw File Name`,
      sample = `Sample ID`,
      metabolite = `Compound Name`,
      area = `Peak Area`
    ) |>
    tidyr::separate(metabolite, c("metabolite", "isotope"), sep = " ") |>
    dplyr::filter(sample != "blank") |>
    dplyr::left_join(samples, by = "sample") |>
    dplyr::filter(stringr::str_detect(end_note, "13c")) |>
    dplyr::group_by(tissue, sample, metabolite) |>
    dplyr::mutate(
      mid = area / sum(area)
    ) |>
    dplyr::select(
      id,
      sample,
      group,
      metabolite,
      isotope,
      area,
      mid
    ) |>
    dplyr::filter(!is.na(mid)) |>
    dplyr::arrange(metabolite, group, sample) |>
    dplyr::group_by(dplyr::across("tissue":"metabolite"))
}
