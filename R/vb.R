# vb.R

read_vb_pharm <- function(x) {
  nms <-
    stringr::str_extract(x, "mct\\d{1}") |>
    toupper()
  purrr::map(
    x,
    \(x) readr::read_csv(
      file = x,
      show_col_types = FALSE
    ) |>
      tidyr::pivot_longer(
        cols = -time,
        names_to = "conc",
        values_to = "value",
        names_transform = as.numeric
      )
  ) |>
    rlang::set_names(nms) |>
    dplyr::bind_rows(.id = "target") |>
    dplyr::mutate(
      target = factor(target, levels = c("MCT4", "MCT1")),
      time = time - 8,
      smooth = slider::slide_dbl(value, mean, .before = 15, .after = 5),
      .by = c(target, conc)
    )
}


plot_vb_rates <- function(z) {
  ggplot2::ggplot(z) +
    ggplot2::facet_wrap(
      ggplot2::vars(target),
      nrow = 1
    ) +
    ggplot2::aes(
      x = time,
      y = value,
      color = conc,
      group = conc
    ) +
    ggplot2::geom_point(
      size = 0.05
    ) +
    ggplot2::scale_color_viridis_c(
      option = "H",
      transform = scales::transform_log(),
      breaks = scales::breaks_log(),
      labels = scales::label_log(),
      oob = scales::oob_squish_infinite
    ) +
    ggplot2::labs(
      x = "Time (s)",
      y = "Fluorescence ratio",
      color = "VB253 (nM)"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.width = ggplot2::unit(0.25, "in"),
    ) +
    NULL
}


read_vb_drc <- function(x) {
  nms <-
    stringr::str_extract(x, "mct\\d{1}") |>
    toupper()
  purrr::map(
    x,
    \(x) readr::read_csv(
      file = x,
      show_col_types = FALSE
    )
  ) |>
    rlang::set_names(nms) |>
    dplyr::bind_rows(.id = "target") |>
    dplyr::group_by(target) |>
    dplyr::mutate(
      target = factor(target, levels = c("MCT4", "MCT1")),
      conc = 10 ^ conc,
      norm = scales::rescale(
        K,
        to = c(0, 1),
        from = c(0, max(K))
      )
    )
}


plot_vb_drc <- function(z) {
  ggplot2::ggplot(z) +
    ggplot2::aes(
      x = conc,
      y = norm,
      color = target,
      fill = target,
      group = target
    ) +
    ggplot2::geom_smooth(
      method = drc::drm,
      method.args = list(fct = drc::L.4()),
      linewidth = 0.5,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      color = "white",
      pch = 21
    ) +
    ggplot2::annotate(
      "text",
      x = 3e-5,
      y = 0.5,
      label = "0.82 nM",
      color = clrs[["MCT4"]],
      size = 8,
      size.unit = "pt",
      family = "Calibri"
    ) +
    ggplot2::annotate(
      "text",
      x = 2,
      y = 0.5,
      label = "26 μM",
      color = clrs[["MCT1"]],
      size = 8,
      size.unit = "pt",
      family = "Calibri"
    ) +
    ggplot2::scale_x_log10(
      breaks = scales::breaks_log(n = 7),
      labels = scales::label_log()
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::labs(
      x = "VB253 (μM)",
      y = "Rate (normalized)",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      legend.position = "bottom"
    )
}


read_vb_invitro_sma <- function(file) {
  readr::read_csv(file, show_col_types = FALSE) |>
    tidyr::pivot_longer(
      tidyselect::matches("_\\d{1}$"),
      names_to = c("measurement", "replicate"),
      names_pattern = "(.+)_(.)",
      values_to = "value"
    ) |>
    tidyr::pivot_wider(
      names_from = "measurement",
      values_from = "value"
    ) |>
    dplyr::mutate(
      type = ifelse(stringr::str_detect(donor, "IPF"), "IPF", "Ctl"),
      .before = donor
    ) |>
    dplyr::group_by(drug, donor) |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(c("sma", "nuclei")),
        \(x) x / mean(x[dose == min(dose)]),
        .names = "{.col}_norm"
      )
    ) |>
    dplyr::group_by(drug, dose, type, donor) |>
    dplyr::summarise(dplyr::across(
      tidyselect::any_of(c("sma", "sma_norm", "nuclei", "nuclei_norm")), mean
    )) |>
    dplyr::mutate(
      drug = factor(
        drug,
        levels = c("VB253", "Nintedanib"),
        labels = c("VB253", "NIN")
      )
    )
}

plot_vb_invitro_sma <- function(
    df,
    source = "IPF",
    y = "sma_norm",
    ytitle = "α-SMA fold change"
) {
  df |>
    dplyr::filter(type == source) |>
    dplyr::group_by(dose) |>
    dplyr::mutate(width = 0.025 * dplyr::n()) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = dose,
      y = .data[[y]],
      color = drug
    ) +
    ggplot2::geom_smooth(
      method = drc::drm,
      formula = y ~ x,
      method.args = list(
        fct = drc::L.4()
      ),
      linewidth = 0.5,
      se = FALSE,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        color = .data$drug,
        width = .data$width
      ),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      # width = 0.1,
      linewidth = 0.25,
      position = ggplot2::position_dodge(width = 0.1),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$drug),
      stat = "summary",
      fun = "mean",
      pch = 21,
      size = 1.5,
      color = "white",
      position = ggplot2::position_dodge(width = 0.1),
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = "Dose (μM)",
      y = ytitle
    ) +
    ggplot2::scale_x_continuous(
      labels = scales::label_math()
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 6, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot()
}

read_vb_mice <- function(file) {
  readr::read_csv(file, show_col_types = FALSE) |>
    refactor() |>
    tidyr::pivot_longer(
      c("ashcroft", "sma", "Penh", "lactate"),
      names_to = "measurement",
      values_to = "value"
    ) |>
    dplyr::filter(!is.na(value))
}
