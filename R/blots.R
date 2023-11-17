# blots.R

clean_blots <- function(df) {
  new_names <-
    c(
      "sma" = "α-SMA",
      "pampk" = "pAMPK",
      "ampk" = "AMPK",
      "psmad3" = "pSmad3",
      "smad3" = "Smad3",
      "erk" = "ERK",
      "perk" = "pERK",
      "lactylk" = "Kla",
      "mct1" = "MCT1",
      "mct4" = "MCT4",
      "col1a1" = "Col1a1",
      "fn1" = "FN1",
      "col3a1" = "Col3a1",
      "cnn1" = "CNN1",
      "hif1a" = "HIF-1α"
    )

  df |>
    dplyr::group_by(dplyr::across(c("experiment":"trial", "protein"))) |>
    dplyr::mutate(
      dplyr::across(
        "background":"density",
        \(x) x / mean(x, na.rm = TRUE)
      ),
      ratio = .data$density / .data$background
    ) |>
    dplyr::group_by(dplyr::across("experiment":"trial")) |>
    dplyr::mutate(protein = new_names[.data$protein]) |>
    refactor() |>
    dplyr::mutate(group = dplyr::cur_group_id()) |>
    dplyr::ungroup()
}
