# contraction.R

# contraction.R

clean_contraction <- function(df) {
  experiment <- batch <- trial <- exp_group <- NULL

  df |>
    dplyr::group_by(.data$experiment, .data$batch, .data$trial) |>
    dplyr::mutate(exp_group = dplyr::cur_group_id(), .after = "trial") |>
    tidyr::pivot_longer(
      tidyselect::matches("d\\d{1}"),
      names_to = "time",
      values_to = "area"
    ) |>
    dplyr::mutate(
      time = stringr::str_extract(.data$time, "\\d"),
      time = as.numeric(.data$time) * 24
    ) |>
    dplyr::group_by(dplyr::across(c("experiment":"treatment", "time"))) |>
    wmo::remove_nested_outliers("area", remove = TRUE) |>
    dplyr::ungroup() |>
    tidyr::complete(
      tidyr::nesting(experiment, batch, trial, exp_group),
      .data$condition, .data$treatment, .data$time
    ) |>
    dplyr::group_by(dplyr::across(c("experiment":"treatment", "time"))) |>
    dplyr::summarise(area = mean(.data$area, na.rm = TRUE)) |>
    dplyr::mutate(frac_start = .data$area / .data$area[.data$time == 0]) |>
    dplyr::filter(.data$time != 0) |>
    dplyr::group_by(dplyr::across(c("experiment":"exp_group", "treatment", "time"))) |>
    dplyr::mutate(frac_ctl = 1 - .data$area / .data$area[.data$condition == "control"]) |>
    dplyr::ungroup() |>
    refactor() |>
    dplyr::mutate(dplyr::across(c("condition", "treatment"), forcats::fct_drop)) |>
    dplyr::rename(group = "exp_group") |>
    dplyr::arrange(.data$experiment, .data$batch, .data$trial, .data$condition, .data$treatment) |>
    dplyr::filter(.data$condition == "TGFβ") |>
    dplyr::filter(.data$time == 24) |>
    dplyr::group_by(.data$treatment) |>
    wmo::remove_nested_outliers("frac_ctl", remove = TRUE) |>
    dplyr::ungroup()
}
