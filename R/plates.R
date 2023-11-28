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

clean_plates_interp <- function(df, cf) {
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

norm_plates <- function(df, cf) {
  dplyr::left_join(
    df,
    cf,
    by = c("experiment", "batch", "trial", "condition", "treatment"),
    suffix = c(".m", ".c")
  ) |>
    dplyr::select(-c(tidyselect::starts_with("assay"), "group.c")) |>
    dplyr::rename(group = "group.m") |>
    dplyr::mutate(norm = .data$conc.m / .data$conc.c)
}

filter_plates <- function(df, exp, treat = NULL, cond = NULL) {
  if (is.null(treat)) {
    treat <- unique(df$treatment)
  }
  if (is.null(cond)) {
    cond <- unique(df$condition)
  }

  df |>
    dplyr::filter(
      .data$experiment %in% exp &
        .data$treatment %in% treat &
        .data$condition %in% cond
    ) |>
    dplyr::mutate(
      dplyr::across(
        c("treatment", "condition"),
        \(x) forcats::fct_drop(x)
      )
    )
}
