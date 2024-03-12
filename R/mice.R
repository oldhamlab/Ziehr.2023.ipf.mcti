# mice.R

# weight ------------------------------------------------------------------

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
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(group = .data$group),
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
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
    ggplot2::geom_text(
      data = stat,
      ggplot2::aes(
        label = label,
        x = 14.5,
        y = y,
        vjust = vjust
      ),
      size = 10 / ggplot2::.pt,
      family = "Calibri"
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

# vent --------------------------------------------------------------------

format_vent <- function(df) {
  df |>
    dplyr::filter(death == 0) |>
    dplyr::select(
      idx,
      start_date,
      death_day,
      genotype,
      source,
      age,
      sex,
      group,
      tidyselect::contains("el_prime"),
      cst
    ) |>
    tidyr::pivot_longer(
      c(tidyselect::contains("el_prime"), cst),
      names_to = c("measurement"),
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(
      replicate = stringr::str_extract(measurement, "\\d{1}$"),
      measurement = stringr::str_extract(measurement, "prime3|prime8|cst"),
      group = factor(
        group,
        levels = c("Veh", "Bleo", "AZD", "VB"),
        labels = c("Ctl", "Bleo", "AZD", "VB"),
        ordered = TRUE
      ),
      vent_date = start_date + lubridate::days(death_day)
    ) |>
    dplyr::group_by(dplyr::across(c(idx:measurement, vent_date))) |>
    dplyr::summarise(value = mean(value)) |>
    dplyr::relocate(vent_date, .after = start_date) |>
    dplyr::filter(idx %nin% c(75, 84)) |>
    dplyr::group_by(.data$measurement) |>
    normalize("value", "start_date") |>
    dplyr::ungroup()
}

stats_vent <- function(df) {
  df |>
    dplyr::group_by(measurement) |>
    tidyr::nest() |>
    dplyr::mutate(stats = purrr::map(data, \(x) {
      lmerTest::lmer(value ~ group + (1 | vent_date), data = x) |>
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
    })) |>
    dplyr::select(measurement, stats) |>
    tidyr::unnest(c(stats))
}

# histo -------------------------------------------------------------------

format_mice_col <- function(df, column = "Ashcroft_avg") {
  x <- rlang::sym(column)
  df <- dplyr::filter(df, !is.na(!!x))

  if (column == "OHP_total") {
    df <- dplyr::filter(df, start_date >= "2022-04-01")
    measure <- "ohp"
  } else if (column == "Ashcroft_avg") {
    measure <- "ashcroft"
  }

  df |>
    dplyr::select(
      idx,
      start_date,
      genotype,
      source,
      age,
      sex,
      group,
      value = x
    ) |>
    dplyr::mutate(
      group = factor(
        group,
        levels = c("Veh", "Bleo", "AZD", "VB"),
        labels = c("Ctl", "Bleo", "AZD", "VB")
      ),
      measurement = measure
    )
}

stats_histo <- function(df, measure, mixed = FALSE, comps = comps_bleo_2) {
  if (mixed) {
    m <- lmerTest::lmer(value ~ group + (1 | animal), data = df)
  } else {
    m <- stats::lm(value ~ group, data = df)
  }

  m |>
    emmeans::emmeans(~ group) |>
    emmeans::contrast(
      method = comps,
      adjust = "dunnettx"
    ) |>
    tibble::as_tibble() |>
    dplyr::rename(group = "contrast") |>
    dplyr::mutate(measurement = measure) |>
    annot_stats() |>
    refactor()
}



# plots -------------------------------------------------------------------

plot_mice <- function(
    df,
    stats,
    measure,
    x = "group",
    y = "value",
    ytitle = NULL
) {
  p <-
    df |>
    dplyr::filter(measurement == measure) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]],
      color = .data[[x]],
      fill = .data[[x]]
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(color = .data$group),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.1,
      linewidth = 0.25,
      position = ggplot2::position_nudge(x = 0.25),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$group),
      stat = "summary",
      fun = "mean",
      pch = 21,
      size = 1.5,
      color = "white",
      position = ggplot2::position_nudge(x = 0.25),
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = ytitle
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

  if (!is.null(stats)) {
    p <-
      p +
      ggplot2::geom_text(
        data = dplyr::filter(stats, measurement == measure),
        ggplot2::aes(
          y = y,
          vjust = vjust,
          label = label
        ),
        color = "black",
        size = 10 / ggplot2::.pt,
        family = "Calibri"
      )
  }
  p
}
