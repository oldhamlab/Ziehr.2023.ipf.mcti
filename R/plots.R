# plots.R

theme_plot <- function(
    base_size = 10,
    base_family = "Arial",
    base_line_size = base_size / 40,
    base_rect_size = base_size / 40
) {

  # The half-line (base-fontsize / 2) sets up the basic vertical
  # rhythm of the theme. Most margins will be set to this value.
  # However, when we work with relative sizes, we may want to multiply
  # `half_line` with the appropriate relative size. This applies in
  # particular for axis tick sizes. And also, for axis ticks and
  # axis titles, `half_size` is too large a distance, and we use `half_size/2`
  # instead.
  half_line <- base_size / 2

  # Throughout the theme, we use three font sizes, `base_size` (`rel(1)`)
  # for normal, `rel(0.8)` for small, and `rel(1.2)` for large.

  theme(
    # elements inherited by others
    line = element_line(
      colour = "grey30",
      linewidth = base_line_size,
      linetype = 1,
      lineend = "butt"
    ),
    rect = element_rect(
      fill = "white",
      colour = NA,
      linewidth = base_rect_size,
      linetype = 1
    ),
    text = element_text(
      family = base_family,
      face = "plain",
      colour = "black",
      size = base_size,
      lineheight = 0.9,
      hjust = 0.5,
      vjust = 0.5,
      angle = 0,
      margin = margin(),
      debug = FALSE
    ),

    axis.line = element_blank(),
    axis.line.x = NULL,
    axis.line.y = NULL,

    axis.text = element_text(
      size = rel(0.8),
      colour = "grey30"
    ),
    axis.text.x = element_text(
      margin = margin(t = 0.8 * half_line / 2),
      vjust = 1
    ),
    axis.text.x.top = element_text(
      margin = margin(b = 0.8 * half_line / 2),
      vjust = 0
    ),
    axis.text.y = element_text(
      margin = margin(r = 0.8 * half_line / 2),
      hjust = 1
    ),
    axis.text.y.right = element_text(
      margin = margin(l = 0.8 * half_line / 2),
      hjust = 0
    ),

    axis.ticks = element_line(colour = "grey30"),
    axis.ticks.length = unit(half_line / 2, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL,

    axis.title.x = element_text(
      margin = margin(t = base_size),
      vjust = 1
    ),
    axis.title.x.top = element_text(
      margin = margin(b = base_size),
      vjust = 0
    ),
    axis.title.y = element_text(
      angle = 90,
      margin = margin(r = base_size),
      vjust = 1
    ),
    axis.title.y.right = element_text(
      angle = -90,
      margin = margin(l = base_size),
      vjust = 0
    ),

    legend.background = element_rect(colour = NA),
    legend.spacing = unit(2 * half_line, "pt"),
    legend.spacing.x = NULL,
    legend.spacing.y = NULL,
    legend.margin = margin(
      half_line,
      half_line,
      half_line,
      half_line
    ),
    legend.key = element_rect(fill = "white", colour = NA),
    legend.key.size = unit(0.8, "lines"),
    legend.key.height = NULL,
    legend.key.width = NULL,
    legend.text = element_text(size = rel(0.8)),
    legend.text.align = NULL,
    legend.title = element_text(hjust = 0),
    legend.title.align = NULL,
    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    legend.box = NULL,
    legend.box.margin = margin(0, 0, 0, 0, "cm"),
    legend.box.background = element_blank(),
    legend.box.spacing = unit(2 * half_line, "pt"),

    panel.background = element_rect(fill = "white", colour = NA),
    panel.border = element_rect(fill = NA, color = "grey30"),
    panel.grid = element_blank(),
    panel.grid.minor = element_line(linewidth = rel(0.5)),
    panel.spacing = unit(half_line, "pt"),
    panel.spacing.x = NULL,
    panel.spacing.y = NULL,
    panel.ontop    = FALSE,

    strip.background = element_rect(fill = "white", colour = NA),
    strip.clip = "inherit",
    strip.text = element_text(
      colour = "grey10",
      size = rel(0.8),
      margin = margin(
        0.8 * half_line,
        0.8 * half_line,
        0.8 * half_line,
        0.8 * half_line
      )
    ),
    strip.text.x = NULL,
    strip.text.y = element_text(angle = -90),
    strip.text.y.left = element_text(angle = 90),
    strip.placement = "inside",
    strip.placement.x = NULL,
    strip.placement.y = NULL,
    strip.switch.pad.grid = unit(half_line / 2, "pt"),
    strip.switch.pad.wrap = unit(half_line / 2, "pt"),

    plot.background = element_rect(colour = "white"),
    plot.title = element_text( # font size "large"
      size = rel(1.2),
      hjust = 0.5,
      vjust = 1,
      margin = margin(b = half_line)
    ),
    plot.title.position = "panel",
    plot.subtitle =      element_text( # font size "regular"
      hjust = 0, vjust = 1,
      margin = margin(b = half_line)
    ),
    plot.caption = element_text( # font size "small"
      size = rel(0.8),
      hjust = 1, vjust = 1,
      margin = margin(t = half_line)
    ),
    plot.caption.position = "panel",
    plot.tag = element_text(
      size = rel(1.2),
      hjust = 0.5, vjust = 0.5
    ),
    plot.tag.position = "topleft",
    plot.margin = margin(half_line, half_line, half_line, half_line),

    complete = TRUE
  )
}


nice_limits <- function(x, n = 5) {
  range(scales::breaks_extended(n = n, only.loose = TRUE)(x))
}


clrs <- c(
  "Ctl"           = "#A0CBE8",
  "IPF"           = "#4E79A7",
  "Bleo"          = "#4E79A7",
  "TGFβ"          = "#4E79A7",
  "Veh"           = "#4E79A7",
  "AZD"           = "#59A14F",
  "VB"            = "#B6992D",
  "AZD/VB"        = "#D37295",
  "Ctl\nNone"     = "#D4A6C8",
  "TGFβ\nNone"    = "#B07AA1",
  "Ctl\nVeh"      = "#A0CBE8",
  "TGFβ\nVeh"     = "#4E79A7",
  "Ctl\nAZD"      = "#8CD17D",
  "TGFβ\nAZD"     = "#59A14F",
  "Ctl\nVB"       = "#F1CE63",
  "TGFβ\nVB"      = "#B6992D",
  "Ctl\nAZD/VB"   = "#FABFD2",
  "TGFβ\nAZD/VB"  = "#D37295",
  NULL
)


plot_one_factor <- function(
    df,
    stats = NULL,
    x,
    y,
    title = NULL,
    ytitle
) {
  p <-
    ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]],
    ) +
    ggplot2::geom_bar(
      stat = "summary",
      fun = "mean",
      ggplot2::aes(fill = .data[[x]]),
      width = 0.6,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.25,
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.75,
      stroke = 0.3,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = ytitle,
      fill = NULL,
      title = title
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 5, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_fill_manual(
      values = clrs,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL

  if (!is.null(stats)) {
    p <-
      p +
      ggplot2::geom_text(
        data = stats,
        ggplot2::aes(
          x = .data[["x"]],
          y = .data[["y"]],
          vjust = .data[["vjust"]],
          label = .data[["label"]]
        ),
        size = 12 / ggplot2::.pt,
        color = "black",
        show.legend = FALSE
      )
  }
  p
}


plot_two_factor <- function(
    df,
    stats = NULL,
    x,
    y,
    title = NULL,
    ytitle,
    wrap = NULL
) {
  dodge_width <- 0.7

  p <-
    ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data[[x]],
      y = .data[[y]],
    ) +
    ggplot2::geom_bar(
      stat = "summary",
      fun = "mean",
      ggplot2::aes(fill = .data[["condition"]]),
      position = ggplot2::position_dodge(width = dodge_width),
      width = 0.6,
      show.legend = TRUE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(group = .data[["condition"]]),
      position = ggplot2::position_dodge(width = dodge_width),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.25,
      linewidth = 0.4,
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(group = .data[["condition"]]),
      dodge.width = dodge_width,
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.75,
      stroke = 0.3,
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = NULL,
      y = ytitle,
      fill = NULL,
      title = title
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::guides(color = "none") +
    ggplot2::scale_fill_manual(
      values = clrs,
      aesthetics = c("fill", "color")
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL

  if (!is.null(stats)) {
    p <-
      p +
      ggplot2::geom_text(
        data = dplyr::filter(stats, is.na(.data[["condition"]])),
        ggplot2::aes(
          y = .data[["y"]],
          label = .data[["label"]],
          vjust = .data[["vjust"]]
        ),
        size = 12 / ggplot2::.pt,
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
        size = 12 / ggplot2::.pt,
        position = ggplot2::position_dodge(width = dodge_width),
        show.legend = FALSE
      )
  }

  if (!is.null(wrap)) {
    p <- p + ggplot2::facet_wrap(ggplot2::vars(.data[[wrap]]))
  }
  p
}
