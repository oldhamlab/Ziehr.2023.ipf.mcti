# ros.R

clean_ros <- function(df) {
  out <-
    df |>
    dplyr::group_by(.data$experiment) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(c("value", "MitoSOX", "MitoTracker")),
        \(x) x - max(0, mean(x[.data$type == "blank"]), na.rm = TRUE)
      )
    ) |>
    dplyr::filter(.data$type != "blank") |>
    dplyr::select(-"type") |>
    refactor() |>
    dplyr::ungroup()

  if ("value" %nin% names(out)) {
    out |>
      dplyr::mutate(
        Ratio = .data$MitoSOX / .data$MitoTracker
      ) |>
      tidyr::pivot_longer(
        c("MitoSOX", "MitoTracker", "Ratio"),
        names_to = "measurement",
        values_to = "value"
      )
  } else if ("value" %in% names(out)) {
    out |>
      dplyr::mutate(measurement = "CellROX", .before = "value")
  }
}
