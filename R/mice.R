# mice.R

# weights -----------------------------------------------------------------

format_wt <- function(df) {
  df |>
    dplyr::filter(.data$death == 0) |>
    dplyr::select(
      idx,
      start_date,
      genotype,
      source,
      age,
      sex,
      group,
      tidyselect::starts_with("wt.")
    ) |>
    tidyr::pivot_longer(
      tidyselect::starts_with("wt."),
      names_to = "day",
      values_to = "wt",
      names_prefix = "wt\\.",
      names_transform = list(day = as.integer)
    ) |>
    dplyr::filter(!is.na(wt)) |>
    dplyr::group_by(idx) |>
    dplyr::mutate(
      wt_norm = wt / wt[day == 0],
      group = factor(
        group,
        levels = c("Veh", "Bleo", "AZD", "VB"),
        labels = c("Ctl", "Bleo", "AZD", "VB"),
        ordered = TRUE
      )
    ) |>
    dplyr::ungroup()
}

bin_wt <- function(df) {
  df |>
    dplyr::mutate(
      wk = round(day / 7) * 7,
      diff = abs(day - wk),
      .after = day
    ) |>
    dplyr::filter(diff == min(diff), .by = c("idx", "wk")) |>
    dplyr::select(-c("day", "diff")) |>
    dplyr::rename(day = "wk")
}

stat_bin_wt <- function(df) {
  lmerTest::lmer(
    wt_norm ~ factor(day, ordered = TRUE) * group + (1 | idx),
    data = df
  ) |>
    emmeans::emmeans(~ group) |>
    emmeans::contrast(
      method = list(
        Ctl = c(-1, 1, 0, 0),
        AZD = c(0, -1, 1, 0),
        VB =  c(0, -1, 0, 1)
      ),
      adjust = "dunnettx"
    ) |>
    tibble::as_tibble() |>
    dplyr::rename(group = contrast) |>
    refactor() |>
    annot_stats()

}

plot_wt_timeline <- function(df, stat) {
  df |>
    dplyr::filter(day != 0) |>
    ggplot2::ggplot() +
    ggplot2::facet_wrap(ggplot2::vars(group), nrow = 1) +
    ggplot2::aes(
      x = day,
      y = wt_norm
    ) +
    ggplot2::geom_hline(
      yintercept = 1,
      linewidth = 0.1,
      linetype = 2
      ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(color = .data$group),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.1,
      linewidth = 0.25,
      position = ggplot2::position_nudge(x = 1),
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      ggplot2::aes(color = .data$group),
      stat = "summary",
      fun = "mean",
      position = ggplot2::position_nudge(x = 1),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$group),
      stat = "summary",
      fun = "mean",
      pch = 21,
      size = 1.5,
      color = "white",
      position = ggplot2::position_nudge(x = 1),
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(group = .data$group),
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = stat,
      ggplot2::aes(
        label = label,
        x = 14.5,
        y = y,
        vjust = vjust
      )
    ) +
    ggplot2::labs(
      x = "Time (d)",
      y = "Weight (fraction of baseline)",
      fill = NULL,
      color = NULL
    ) +
    ggplot2::scale_x_continuous(
      breaks = seq(7, 21, 7),
      limits = c(7, 22)
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}
