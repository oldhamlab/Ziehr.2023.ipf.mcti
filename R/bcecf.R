# bcecf.R

format_bcecf <- function(x) {
  x |>
    tidyr::pivot_longer(
      `1`:tidyselect::last_col(),
      names_to = "col",
      values_to = "value",
      names_transform = list(col = as.numeric)
    ) |>
    dplyr::mutate(
      type = dplyr::case_when(
        col %in% 1:2 ~ "blank",
        col %in% 3:10 ~ "sample",
        col %in% 11:12 ~ "standard"
      ),
      condition = dplyr::case_when(
        col %in% 3:6 ~ "Ctl",
        col %in% 7:10 ~ "TGFβ"
      ),
      condition = factor(condition, levels = c("Ctl", "TGFβ")),
      treatment = dplyr::case_when(
        col %in% c(3, 7) ~ "Veh",
        col %in% c(4, 8) ~ "AZD",
        col %in% c(5, 9) ~ "VB",
        col %in% c(6, 10) ~ "AZD/VB"
      ),
      treatment = factor(treatment, levels = c("Veh", "AZD", "VB", "AZD/VB")),
      group = dplyr::case_when(
        col == 1 ~ "blank",
        col == 2 ~ "cells",
        col %in% 3:10 ~ stringr::str_c(condition, ".", treatment),
        col == 11 & row %in% 1:4 ~ "4.5",
        col == 11 & row %in% 5:8 ~ "5.5",
        col == 12 & row %in% 1:4 ~ "6.5",
        col == 12 & row %in% 5:8 ~ "7.5"
      ),
      .before = value
    )
}

clean_bcecf <- function(x) {
  x |>
    dplyr::filter(col %in% 2:10) |>
    dplyr::group_by(date, excitation) |>
    dplyr::mutate(value = value - mean(value[group == "cells"])) |>
    dplyr::filter(type != "blank") |>
    dplyr::group_by(date, row, col, condition, treatment, group) |>
    dplyr::summarize(ratio = value[excitation == 490] / value[excitation == 440]) |>
    dplyr::group_by(date, condition, treatment, group) |>
    wmo::remove_nested_outliers(ratio, remove = TRUE) |>
    dplyr::summarise(ratio = mean(ratio)) |>
    dplyr::select(group = date, condition, treatment, ratio) |>
    dplyr::ungroup()
}

plot_bcecf <- function(df, stats) {
  normalize(df, "ratio", "group") |>
    dplyr::rename(experiment = "group") |>
    dplyr::mutate(group = stringr::str_c(condition, treatment, sep = "\n")) |>
    refactor() |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = treatment,
      y = ratio_corr,
      color = condition,
      fill = condition
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(color = .data$condition),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.1,
      linewidth = 0.25,
      position = ggplot2::position_dodge(width = 0.8),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$condition),
      stat = "summary",
      fun = "mean",
      pch = 21,
      size = 1.5,
      color = "white",
      position = ggplot2::position_dodge(width = 0.8),
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      dodge.width = 0.8,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(stats, is.na(.data[["condition"]])),
      ggplot2::aes(
        y = .data[["y"]],
        label = .data[["label"]],
        vjust = .data[["vjust"]]
      ),
      size = 10 / ggplot2::.pt,
      family = "Calibri",
      color = "black",
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = dplyr::filter(stats, !is.na(.data[["condition"]])),
      ggplot2::aes(
        y = .data[["y"]],
        label = .data[["label"]],
        vjust = .data[["vjust"]],
        color = .data[["condition"]]
      ),
      size = 10 / ggplot2::.pt,
      family = "Calibri",
      position = ggplot2::position_dodge(width = 0.8),
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Fluorescence ratio"
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 6, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}
