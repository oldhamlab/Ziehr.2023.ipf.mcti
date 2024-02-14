# ros.R

clean_ros <- function(df) {
  x <-
    df |>
    dplyr::group_by(.data$experiment) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(c("value", "MitoSOX", "MitoTracker")),
        \(x) x - max(0, mean(x[.data$type == "blank"]), na.rm = TRUE)
      ),
      group = stringr::str_c(.data$condition, .data$treatment, sep = "\n")
    ) |>
    dplyr::filter(.data$type != "blank") |>
    dplyr::select(-"type") |>
    dplyr::relocate(.data$group, .after = .data$treatment) |>
    refactor() |>
    dplyr::ungroup()

  if ("value" %nin% names(x)) {
    out <-
      x |>
      dplyr::mutate(
        Ratio = .data$MitoSOX / .data$MitoTracker
      ) |>
      tidyr::pivot_longer(
        c("MitoSOX", "MitoTracker", "Ratio"),
        names_to = "measurement",
        values_to = "value"
      )
  } else if ("value" %in% names(x)) {
    out <-
      x |>
      dplyr::mutate(measurement = "CellROX", .before = "value")
  }
  out |>
    dplyr::group_by(
      dplyr::across(c("experiment", "condition":"measurement"))
    ) |>
    wmo::remove_nested_outliers("value", remove = TRUE) |>
    dplyr::summarise(value = mean(.data$value)) |>
    dplyr::group_by(.data$experiment, .data$measurement, .data$time) |>
    dplyr::mutate(value_norm = .data$value / mean(.data$value)) |>
    dplyr::ungroup()
}

norm_ros <- function(df, measure) {
  df |>
    dplyr::filter(measurement == measure) |>
    dplyr::filter(time == 3) |>
    dplyr::mutate(
      value_corr = .data$value_norm / mean(.data$value_norm[.data$group == "TGFβ\nVeh"])
    )
}
