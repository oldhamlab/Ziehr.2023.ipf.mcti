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

make_fig03 <- function(p1, p2, p3, p4, p5) {
  design <- "abc \n ddd \n eff"

  ((p1 + p2 + patchwork::guide_area() & ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))) +
      patchwork::plot_layout(guides = "collect")) +
    ((p3) + patchwork::plot_layout(guides = "keep")) +
    (p4 + patchwork::plot_layout(guides = "keep")) + p5 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(2, 2, 0.5), "in"),
      heights = ggplot2::unit(c(1.5), "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_fig03s <- function(p1, p2, p3, p4) {
  design <- "abc \n ddd \n eee"

  ((p1 + p2 + patchwork::guide_area() & ggplot2::guides(fill = ggplot2::guide_legend(nrow = 2))) +
      patchwork::plot_layout(guides = "collect")) +
    ((p3) + patchwork::plot_layout(guides = "keep")) +
    p4 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(2, 2, 0.5), "in"),
      heights = ggplot2::unit(c(1.5, 2, 1.5), "in")
    )
}

make_fig04 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  design <- "abcd"
  extra <-
    p1 + p2 + p3 + p4 +
    patchwork::plot_annotation(
      title = "Extracellular",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = ggplot2::rel(1)
        )
      )
    ) +
    theme_patchwork(
      design = design,
      tags = list(c("A", "B", "C", "D")),
      widths = ggplot2::unit(2, "in"),
      heights = ggplot2::unit(2, "in")
    )

  intra <-
    p5 + p6 + p7 + p8 +
    patchwork::plot_annotation(
      title = "Intracellular",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = ggplot2::rel(1)
        )
      )
    ) +
    theme_patchwork(
      design = design,
      tags = list(c("E", "F", "G", "H")),
      widths = ggplot2::unit(2, "in"),
      heights = ggplot2::unit(2, "in")
    )

  mids <-
    p9 + p10 +
    patchwork::plot_layout(guides = "collect") +
    theme_patchwork(
      design = "ab",
      tags = list(c("I", "J"))
    ) &
    ggplot2::theme(legend.position = "bottom")

  patchwork::wrap_elements(full = extra) +
    patchwork::wrap_elements(full = intra) +
    patchwork::wrap_elements(full = mids) +
    theme_patchwork(
      design = "a \n b \n c",
      widths = ggplot2::unit(9.5, "in"),
      heights = ggplot2::unit(3.5, "in"),
      tags = NULL
    )
}

make_fig04s <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13) {
  design <- "ab \n cd \n ef"

  extra <-
    p1 + p2 + p5 + p6 + p9 + p10 +
    patchwork::plot_annotation(
      title = "Extracellular",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = ggplot2::rel(0.8),
          face = "bold"
        )
      )
    ) +
    theme_patchwork(
      design = design,
      tags = list(c("A", "B", "E", "F", "I", "J"))
    )
  intra <-
    p3 + p4 + p7 + p8 + p11 + p12 +
    patchwork::plot_annotation(
      title = "Intracellular",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = ggplot2::rel(0.8),
          face = "bold"
        )
      )
    ) +
    theme_patchwork(
      design = design,
      tags = list(c("C", "D", "G", "H", "K", "L"))
    )
  mids <-
    p13 +
    theme_patchwork(
      tags = list("M")
    ) &
    ggplot2::theme(legend.position = "bottom")

  patchwork::wrap_elements(full = extra) +
    patchwork::wrap_elements(full = intra) +
    patchwork::wrap_elements(full = mids) +
    theme_patchwork(
      design = "ab \n c#",
      widths = ggplot2::unit(4, "in"),
      heights = ggplot2::unit(c(6, 3), "in"),
      tags = NULL
    )
}

make_fig05 <- function(p1, p2, p3, p4) {
  design <- "ab \n cb \n dd"
  p1 + p2 + p3 + p4 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(2, 1.5), "in"),
      heights = ggplot2::unit(c(2, 2, 2), "in")
    )
}

make_fig05s <- function(p1, p2, p3) {
  design <- "abc"
  p1 + p2 + p3 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(2, "in"),
      heights = ggplot2::unit(2, "in")
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_fig06 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10) {
  design <- "abcd \n efgh \n ij##"
  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(c(1, 2, 1), "in"),
      guides = "collect"
    ) &
    ggplot2::theme(legend.position = "bottom")
}

make_fig07 <- function(p1, p2, p3, p4) {
  x <- p1 + p2 + p3 + patchwork::plot_layout(guides = "collect")
  x - p4 +
    theme_patchwork(
      design = "aaab",
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(1.5, "in"),
      tags = list(c("A", "", "", "B"))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -10)
    )
}

make_fig07s <- function(p1, p2, p3, p4, p5, p6, p7, p8) {
  p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 +
    theme_patchwork(
      design = "abc \n def \n gh#",
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(1.5, "in"),
      guides = "collect",
      tags = list(c("A", "", "", "B", "", "", "C", "D"))
    ) &
    ggplot2::theme(legend.position = "bottom")
}
