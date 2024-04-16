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
      theme = ggplot2::theme(plot.margin = ggplot2::margin(rep(3, 4)))
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
    units = "in",
    res = 220
  )

  if (file.exists("Rplots.pdf")) unlink("Rplots.pdf")

  stringr::str_c(path, "/", filename, ".png")
}

make_fig01 <- function(p1, p2, p3, p4, p5, p6) {
  figs <- plot_scale(p1, p3, p5)
  design <- "ab \n cd \n ef"
  figs[[1]] + p2 + figs[[2]] + p4 + figs[[3]] + p6 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1.5, 2), "in"),
      heights = ggplot2::unit(1.5, "in"),
      tags = list(c("A", "", "B", "", "C", ""))
    )
}

make_fig02 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) {
  figs <- plot_scale(p1, p3, p5)

  col1 <-
    figs[[1]] + p2 + figs[[2]] + p4 + figs[[3]] + p6 + p7 + p8 +
    theme_patchwork(
      design = "ab \n cd \n ef \n gh",
      widths = c(1.75, 1.5),
      heights = c(2, 2, 2, 1),
      guides = "collect",
      tags = list(c("A", "", "B", "", "C", "", "D", ""))
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )

  col2 <-
    p9 + p10 + p11 + p12 +
    theme_patchwork(
      design = "ab \n cb \n dd",
      widths = c(2, 1.5),
      heights = c(2, 2, 2),
      tags = list(c("E", "F", "G", "H"))
    )

  patchwork::wrap_elements(full = col1) +
    patchwork::wrap_elements(full = col2) +
    theme_patchwork(
      design = "ab",
      widths = ggplot2::unit(c(4, 6), "in"),
      heights = ggplot2::unit(8, "in"),
      tags = NULL
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

make_fig04 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16) {
  figs <- plot_scale(p4, p8)

  design <- "abcd"
  extra <-
    p1 + p2 + p3 + figs[[1]] +
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
      widths = c(1, 1, 1, 1)
    )

  intra <-
    p5 + p6 + p7 + figs[[2]] +
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
      widths = c(1, 1, 1, 1)
    )

  mids <-
    p9 + p10 + p11 +
    patchwork::plot_layout(guides = "collect") +
    theme_patchwork(
      design = "abc",
      widths = c(0.4, 1, 1),
      tags = list(c("I", "J", "K"))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -5)
    )

  ros <-
    ((p12 + p13 + p14) + patchwork::plot_layout(guides = "collect")) +
    p15 + p16 +
    theme_patchwork(
      design = "abcde",
      tags = list(c("L", "", "", "M", "N"))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -10)
    )

  patchwork::wrap_elements(full = extra) +
    patchwork::wrap_elements(full = intra) +
    patchwork::wrap_elements(full = mids) +
    patchwork::wrap_elements(full = ros) +
    theme_patchwork(
      design = "a \n b \n c \n d",
      widths = ggplot2::unit(9.5, "in"),
      heights = ggplot2::unit(c(3, 3, 3, 2), "in"),
      tags = NULL
    )
}

make_fig05 <- function(p1, p2, p3, p4) {
  figs <- plot_scale(p1, p3)
  design <- "ac \n bd"
  figs[[1]] + p2 + figs[[2]] + p4 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(c(1, 2), "in"),
      guides = "collect",
      tags = list(c("A", "", "B", ""))
    ) &
    ggplot2::theme(legend.position = "bottom")
}

make_fig06 <- function(p1, p2, p3, p4, p5, p6) {
  design <- "abc \n def"
  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(2, 1.5, 1.5), "in"),
      heights = ggplot2::unit(1.5, "in")
    )
}

make_fig07 <- function(p1, p2, p3) {
  design <- "aaa \n bbb \n cc#"
  p1 + p2 + p3 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1.25), "in"),
      heights = ggplot2::unit(c(1, 3, 1), "in")
    )
}

make_fig08 <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11) {
  # design <- "abcc \n defg \n hijk"
  # p1 + p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9 + p10 + p11 +
  #   theme_patchwork(
  #     design = design,
  #     widths = ggplot2::unit(1.5, "in"),
  #     heights = ggplot2::unit(1.5, "in")
  #   )

  top <-
    ((p1 + p2) +
       patchwork::plot_layout(guides = "collect") &
       ggplot2::theme(
         legend.position = "bottom",
         legend.box.margin = ggplot2::margin(t = -10)
       )) -
    p3 +
    theme_patchwork(
      # design = "abc",
      widths = c(2, 2),
      tags = list(c("A", "B", "C")),
      guides = "collect"
    ) &
    ggplot2::theme(legend.position = "bottom")

  mid <-
    p4 + p5 + p6 + p7 +
    patchwork::plot_annotation(
      title = "Young mice",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = ggplot2::rel(0.8),
          face = "bold"
        )
      )
    ) +
    theme_patchwork(
      design = "abcd",
      tags = list(c("D", "E", "F", ""))
    )

  bot <-
    p8 + p9 + p10 + p11 +
    patchwork::plot_annotation(
      title = "Aged mice",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          size = ggplot2::rel(0.8),
          face = "bold"
        )
      )
    ) +
    theme_patchwork(
      design = "abcd",
      tags = list(c("G", "H", "", "I"))
    )

  patchwork::wrap_elements(full = top) +
    patchwork::wrap_elements(full = mid) +
    patchwork::wrap_elements(full = bot) +
    theme_patchwork(
      design = "a \n b \n c",
      widths = ggplot2::unit(7.5, "in"),
      heights = ggplot2::unit(2, "in"),
      tags = NULL
    )
}

make_fig01s <- function(p1, p2, p3, p4, p5) {
  design <- "a# \n bc \n de"
  p1 + p2 + p3 + p4 + p5 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(c(1.5, 1.5), "in"),
      heights = ggplot2::unit(1.5, "in"),
      guides = "collect",
      tags = list(c("A", "B", "", "C", "D"))
    ) &
    ggplot2::theme(
      legend.position = "bottom"
    )
}

make_fig02s <- function(p1, p2, p3) {
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



make_fig04s <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13) {
  figs <- plot_scale(p2, p4, p6, p8, p10, p12)

  design <- "ab \n cd \n ef"

  extra <-
    p1 + figs[[1]] + p5 + figs[[3]] + p9 + figs[[5]] +
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
      widths = c(1, 0.8, 1, 0.8),
      tags = list(c("A", "B", "E", "F", "I", "J"))
    )
  intra <-
    p3 + figs[[2]] + p7 + figs[[4]] + p11 + figs[[6]] +
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
      widths = c(1, 0.8, 1, 0.8),
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
      heights = ggplot2::unit(c(7.5, 3), "in"),
      tags = NULL
    )
}

make_fig06s <- function(p1, p2, p3, p4, p5, p6) {
  p1 + p2 + p3 + p4 + p5 + p6 +
    theme_patchwork(
      design = "abc \n def",
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(1.5, "in"),
      guides = "collect",
      tags = list(c("A", "", "", "B", "C", "D"))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -5)
    )
}

make_fig07s <- function(p1, p2, p3, p4) {
  figs <- plot_scale(p1, p3)
  design <- "ac \n bd"
  figs[[1]] + p2 + figs[[2]] + p4 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(1, "in"),
      guides = "collect",
      tags = list(c("A", "", "B", ""))
    ) &
    ggplot2::theme(
      legend.position = "bottom",
      legend.margin = ggplot2::margin(t = -5)
    )
}

make_fig08s <- function(p1) {
  design <- "a"
  p1 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(4.5, "in"),
      heights = ggplot2::unit(1.25, "in"),
      tags = NULL
    )
}



make_fig09s <- function(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12) {
  figs <- plot_scale(p2, p4, p6, p8, p10, p12)
  # figs <- list(p2, p4, p6, p8, p10, p12)

  design <- "ab \n cd \n ef"

  extra <-
    p1 + figs[[1]] + p5 + figs[[3]] + p9 + figs[[5]] +
    patchwork::plot_annotation(
      title = "Plasma",
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
      # widths = c(1, 0.8, 1, 0.8),
      tags = list(c("A", "B", "E", "F", "I", "J"))
    )
  intra <-
    p3 + figs[[2]] + p7 + figs[[4]] + p11 + figs[[6]] +
    patchwork::plot_annotation(
      title = "Lung",
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
      # widths = c(1, 0.8, 1, 0.8),
      tags = list(c("C", "D", "G", "H", "K", "L"))
    )

  patchwork::wrap_elements(full = extra) +
    patchwork::wrap_elements(full = intra) +
    theme_patchwork(
      design = "ab",
      widths = ggplot2::unit(4, "in"),
      heights = ggplot2::unit(7, "in"),
      tags = NULL
    )
}

# make_fig09s2 <- function(p1, p2) {
#   design <- "aab"
#   p1 + p2 +
#     theme_patchwork(
#       design = design,
#       widths = ggplot2::unit(1.5, "in"),
#       heights = ggplot2::unit(2, "in")
#     )
# }

make_fig10s <- function(p1) {
  design <- "a"
  p1 +
    theme_patchwork(
      design = design,
      widths = ggplot2::unit(1.5, "in"),
      heights = ggplot2::unit(1.5, "in"),
      tags = NULL
    )
}

create_resources <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    flextable::flextable() |>
    flextable::bold(part = "header", bold = TRUE) |>
    flextable::font(fontname = "Calibri", part = "all") |>
    flextable::fontsize(size = 9, part = "all") |>
    flextable::set_table_properties(layout = "autofit")
}
