# plates.R

clean_plates <- function(df) {
  df |>
    tidyr::pivot_longer("a":"c", names_to = "replicate", values_to = "value") |>
    dplyr::filter(!is.na(.data$value)) |>
    dplyr::group_by(dplyr::across(c("experiment":"conc"))) |>
    wmo::remove_nested_outliers("value", remove = TRUE) |>
    dplyr::summarise(value = mean(.data$value)) |>
    refactor() |>
    dplyr::ungroup()
}

clean_data <- function(df, cf) {
  if ("well" %in% names(df)) {
    df <-
      df |>
      dplyr::select(-"assay") |>
      dplyr::group_by(dplyr::across(!c("well", "conc"))) |>
      wmo::remove_nested_outliers(column = "conc", remove = TRUE) |>
      dplyr::summarise(conc = mean(.data$conc)) |>
      dplyr::group_by(dplyr::across("experiment":"trial")) |>
      dplyr::mutate(group = dplyr::cur_group_id()) |>
      dplyr::ungroup() |>
      dplyr::arrange(.data$group)
  }
  dplyr::mutate(df, conc = .data$conc * cf)
}
