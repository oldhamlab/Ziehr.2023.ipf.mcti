# seahorse.R

format_wells <- function(file) {
  readr::read_csv(file, show_col_types = FALSE) |>
    refactor() |>
    dplyr::mutate(group = stringr::str_c(condition, treatment, sep = "\n"))
}

plot_seahorse_timelines <- function(x) {
  df <-
    seahorse::rates(x) |>
    dplyr::mutate(
      experiment = stringr::str_extract(experiment, "\\d{4}-\\d{2}-\\d{2}"),
      name = factor(
        rate,
        levels = c("PER", "OCR"),
        labels = c("Proton efflux rate", "Oxygen consumption rate")
      ),
      value = .data$value * 1000
    ) |>
    refactor()

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$measurement,
      y = .data$value,
      fill = .data$group,
      color = .data$group
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$name),
      scales = "free_y"
    ) +
    ggplot2::stat_summary(
      geom = "line",
      fun = "mean",
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "errorbar",
      fun.data = ggplot2::mean_se,
      width = 0.2,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::stat_summary(
      geom = "point",
      fun = "mean",
      pch = 21,
      color = "white",
      size = 1.25
    ) +
    ggplot2::geom_point(alpha = 0) +
    ggplot2::scale_x_continuous(
      breaks = c(3.5, 6.5, 9.5),
      labels = c("Oligo", "FCCP", "Rot/AMA")
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill"),
      limits = force
    ) +
    ggplot2::labs(
      y = "fmol/min/cell",
      x = NULL,
      color = NULL,
      fill = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot()
}

summarize_seahorse <- function(x) {
  seahorse::summary(x) |>
    dplyr::group_by(.data$rate, .data$stage) |>
    dplyr::mutate(
      value = .data$value * 1000,
      grand_mean = mean(.data$value),
      experiment = stringr::str_extract(experiment, "\\d{4}-\\d{2}-\\d{2}"),
      name = factor(
        rate,
        levels = c("PER", "OCR"),
        labels = c("Proton efflux rate", "Oxygen consumption rate")
      ),
      stage = factor(
        stage,
        levels = c("basal", "oligo", "fccp", "rot/ama"),
        labels = c("Basal", "Oligo", "FCCP", "Rot/AMA")
      )
    ) |>
    refactor() |>
    dplyr::group_by(.data$rate, .data$name, .data$stage, .data$experiment) |>
    dplyr::mutate(
      exp_mean = mean(.data$value),
      cf = .data$grand_mean - .data$exp_mean,
      value_corr = .data$value + .data$cf
    ) |>
    dplyr::group_by(.data$rate, .data$name, .data$stage)
}

seahorse_group_onefactor <- function(x) {
  group_onefactor(
    x,
    comps = comps_azd_3,
    fo = stats::formula("value_corr ~ group + (1 | experiment)"),
    emm = stats::formula("~ group")
  ) |>
    dplyr::mutate(
      x = stringr::str_c(.data$condition, .data$treatment, sep = "\n")
    ) |>
    refactor()
}

plot_seahorse_summary <- function(x, stats) {
  plot_one_factor(x, stats, x = "group", y = "value_corr", ytitle = "fmol/min/cell") +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$stage),
      rows = ggplot2::vars(.data$name),
      scales = "free_y"
    )
}
