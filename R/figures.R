# figures.R

theme_patchwork <- function(
    design = NULL,
    widths = NULL,
    heights = NULL,
    tags = "A",
    ...
) {
  list(
    patchwork::plot_layout(
      design = design,
      widths = widths,
      heights = heights,
      ...
    ),
    patchwork::plot_annotation(
      tag_levels = tags,
      theme = ggplot2::theme(plot.margin = ggplot2::margin(5, 5, 5, 5))
    )
  )
}

write_figures <- function(plot, filename, path = "manuscript/figs") {
  gtab <- patchwork::patchworkGrob(plot)

  overall_width <-
    grid::convertWidth(
      sum(gtab$widths),
      unitTo = "in",
      valueOnly = TRUE
    )

  overall_height <-
    grid::convertHeight(
      sum(gtab$heights),
      unitTo = "in",
      valueOnly = TRUE
    )

  ggplot2::ggsave(
    filename = stringr::str_c(filename, ".png"),
    plot = plot,
    device = ragg::agg_png,
    path = path,
    width = overall_width,
    height = overall_height,
    units = "in"
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename, ".png")
}

make_fig01 <- function(p1, p2, p3, p4, p5, p6) {
  design <- "ab \n cd \n ef"
  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1, 1.75), "in"),
      heights = ggplot2::unit(1, "in")
    )
}

make_fig02 <- function(p1, p2, p3, p4, p5, p6, p7, p8) {
  design <- "
  abcd
  abcd
  efgh
  ef##
  "

  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1.5, 1.75, 1.5, 1.75), "in"),
      heights = ggplot2::unit(c(0.75, 0.75, 1, 0.5), "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_fig02s <- function(p1, p2, p3, p4, p5) {
  design <- "a# \n bc \n de"
  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1.5, 1.5), "in"),
      heights = ggplot2::unit(1, "in"),
      guides = "collect"
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}
