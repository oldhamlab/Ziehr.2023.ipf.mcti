# nad.R

clean_nad <- function(protein, nad) {
  protein <-
    protein |>
    dplyr::mutate(
      protein = 300 * .data$conc  # μg protein
    ) |>
    dplyr::select(-"conc")

  nad <-
    nad |>
    dplyr::rename(measurement = "nucleotide") |>
    dplyr::mutate(
      conc = 900 * .data$conc  # pmol nucleotide
    )

  dplyr::left_join(
    protein,
    nad,
    by = c("experiment", "sample", "condition", "treatment")
  ) |>
    tidyr::pivot_wider(
      names_from = "measurement",
      values_from = "conc"
    ) |>
    dplyr::mutate(
      NAD_norm = .data$NAD / .data$protein,  # pmol / μg
      NADH_norm = .data$NADH / .data$protein,
      Ratio = .data$NADH / .data$NAD
    ) |>
    refactor() |>
    tidyr::pivot_longer(
      "NAD":"Ratio",
      names_to = "measurement",
      values_to = "value"
    ) |>
    dplyr::rename(group = "experiment") |>
    dplyr::ungroup()
}

norm_nad <- function(df) {
  df |>
    dplyr::mutate(grand_mean = mean(.data$value)) |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(
      exp_mean = mean(.data$value),
      cf = .data$grand_mean - .data$exp_mean,
      value_corr = .data$value + .data$cf
    )
}
