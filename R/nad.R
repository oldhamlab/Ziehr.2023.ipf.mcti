# nad.R

clean_nad <- function(protein, nad) {
  metabs <- unique(nad$nucleotide)
  top <- metabs[metabs %in% c("NADH", "NADPH")]
  bottom <- metabs[metabs %in% c("NAD", "NADP")]

  protein <-
    protein |>
    dplyr::rename(protein = "conc") |>
    # μg protein or divide cell by 1000
    dplyr::mutate(protein = ifelse(top == "NADPH", protein / 1000, protein))

  nad <-
    nad |>
    dplyr::rename(measurement = "nucleotide") # pmol nucleotide

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
      # pmol / μg protein or fmol / cell
      "{bottom}_norm" := .data[[bottom]] / .data$protein,
      "{top}_norm" := .data[[top]] / .data$protein,

      "{top}/{bottom}" := .data[[top]] / .data[[bottom]]
    ) |>
    refactor() |>
    tidyr::pivot_longer(
      !("experiment":"protein"),
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
