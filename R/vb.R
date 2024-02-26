# vb.R

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
