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
      )
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
    dplyr::group_by(.data$rate, .data$name, .data$stage) |>
    normalize("value", "experiment")
}

seahorse_group_onefactor <- function(x) {
  group_onefactor(
    x,
    comps = comps_azd_3,
    mixed = FALSE,
    fo = stats::formula("value_corr ~ group"),
    emm = stats::formula("~ group")
  ) |>
    dplyr::mutate(
      x = stringr::str_c(.data$condition, .data$treatment, sep = "\n")
    ) |>
    refactor()
}

plot_seahorse_summary <- function(x, stats) {
  plot_one_factor(
    x,
    stats,
    x = "group",
    y = "value_corr",
    ytitle = "fmol/min/cell"
  ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(.data$stage),
      rows = ggplot2::vars(.data$name),
      scales = "free_y"
    )
}

summarize_seahorse_atp <- function(x) {
  seahorse::atp(x) |>
    tidyr::separate(
      .data$group,
      c("condition", "treatment"),
      sep = "\n",
      remove = FALSE
    ) |>
    dplyr::group_by(.data$rate) |>
    normalize("value", "experiment") |>
    refactor()
}

plot_seahorse_atp_bars <- function(x, stats) {
  df <-
    x |>
    dplyr::select("experiment", "group", "rate", "value_corr")

  means <-
    df |>
    tidyr::pivot_wider(names_from = "rate", values_from = "value_corr") |>
    dplyr::group_by(.data$group) |>
    dplyr::summarize(
      mito = ggplot2::mean_se(.data$mito),
      glyco = ggplot2::mean_se(.data$glyco)
    ) |>
    tidyr::unnest_wider(
      c("mito", "glyco"),
      names_sep = "_"
    ) |>
    dplyr::mutate(
      dplyr::across(tidyselect::contains("glyco_ym"), \(x) .data$mito_y + x)
    ) |>
    tidyr::pivot_longer(
      !"group",
      names_to = c("rate", "measure"),
      names_sep = "_"
    ) |>
    tidyr::pivot_wider(names_from = "measure", values_from = "value") |>
    dplyr::left_join(stats, by = c("group" = "x", "rate"), suffix = c("", ".stat")) |>
    dplyr::group_by(.data$group) |>
    dplyr::mutate(
      y.stat = dplyr::case_when(
        .data$rate == "mito" ~ y,
        .data$rate == "glyco" ~ .data$y + .data$y[.data$rate == "mito"]
      ),
      rate = stringr::str_to_title(.data$rate)
    ) |>
    refactor()

  pts <-
    df |>
    dplyr::mutate(rate = stringr::str_to_title(.data$rate)) |>
    dplyr::left_join(
      means |>
        dplyr::filter(.data$rate == "Mito") |>
        dplyr::select("group", "rate", "y"),
      by = "group"
    ) |>
    dplyr::mutate(
      value = ifelse(
        .data$rate.x == .data$rate.y,
        .data$value_corr,
        .data$value_corr + .data$y
      )
    ) |>
    refactor()

  ggplot2::ggplot(means) +
    ggplot2::aes(
      x = .data$group,
      y = .data$y
    ) +
    ggplot2::geom_col(
      ggplot2::aes(fill = .data$rate),
      width = 0.6
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$ymin,
        ymax = .data$ymax
      ),
      width = 0.2,
      linewidth = 0.3,
    ) +
    ggbeeswarm::geom_quasirandom(
      data = pts,
      ggplot2::aes(
        y = .data$value,
        shape = .data$rate.x
      ),
      color = "black",
      width = 0.05,
      size = 0.75,
      stroke = 0.2,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        y = .data$y.stat,
        label = .data$label
      ),
      size = 10 / ggplot2::.pt,
      family = "Calibri",
      nudge_y = 2
    ) +
    ggplot2::labs(
      x = NULL,
      y = "ATP production (fmol/cell/min)",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 5, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_fill_brewer(
      palette = "Set1",
      aesthetics = c("fill", "color"),
      direction = -1
    ) +
    ggplot2::scale_shape_manual(
      values = c(1, 2)
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      legend.box.margin = ggplot2::margin(t = -30)
    ) +
    NULL
}

plot_seahorse_pheno <- function(x) {
  x |>
    dplyr::select("experiment", "group", "rate", "value_corr") |>
    tidyr::pivot_wider(names_from = "rate", values_from = "value_corr") |>
    dplyr::group_by(.data$group) |>
    dplyr::summarise(
      mito = ggplot2::mean_se(.data$mito),
      glyco = ggplot2::mean_se(.data$glyco)
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data$mito$y,
      y = .data$glyco$y,
      xmin = .data$mito$ymin,
      xmax = .data$mito$ymax,
      ymin = .data$glyco$ymin,
      ymax = .data$glyco$ymax,
      fill = .data$group,
      color = .data$group
    ) +
    ggplot2::geom_errorbar(
      width = 0.25,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbarh(
      height = 0.25,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 1.5,
      pch = 21,
      color = "white",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      ggplot2::aes(
        label = .data$group
      ),
      point.padding = 6,
      force_pull = 1,
      lineheight = 1,
      nudge_x = 0.5,
      nudge_y = -0.5,
      color = "black",
      family = "Calibri",
      size = 1.75,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_extended(only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::labs(
      x = "Mito ATP",
      y = "Glyco ATP",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::coord_fixed(
      expand = TRUE,
      clip = "off"
    ) +
    theme_plot() +
    NULL
}

summarize_seahorse_stress <- function(x) {
  seahorse::mst(x) |>
    dplyr::bind_rows(seahorse::gst(x)) |>
    tidyr::separate(
      .data$group,
      c("condition", "treatment"),
      sep = "\n",
      remove = FALSE
    ) |>
    dplyr::group_by(.data$rate) |>
    dplyr::mutate(
      rate = factor(
        rate,
        levels = c("gst", "src", "coupling"),
        labels = c(
          "Glycolytic capacity",
          "Spare respiratory capacity",
          "Coupling efficiency"
        )
      )
    ) |>
    normalize("value", "experiment") |>
    refactor()
}

plot_seahorse_stress <- function(x, stats) {
  plot_one_factor(x, stats, x = "group", y = "value_corr", ytitle = "Ratio") +
    ggplot2::facet_wrap(
      facets = ggplot2::vars(.data$rate),
      nrow = 1,
      scales = "free_y"
    )
}
