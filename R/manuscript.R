# manuscript.R

write_pkg_cites <- function() {
  desc::desc_get_deps() |>
    dplyr::pull(.data$package) |>
    knitr::write_bib(file = manuscript_path("pkgs.bib")) |>
    suppressWarnings()
}

test_quarto <- function(file = "manuscript.qmd") {
  quarto::quarto_render(
    manuscript_path(file),
    execute_dir = getwd()
  )
}

plot_image <- function(img, scale = 1, hjust = 0, vjust = 0) {
  cowplot::ggdraw() +
    cowplot::draw_image(
      magick::image_read(img),
      scale = scale,
      hjust = hjust,
      vjust = vjust
    )  +
    theme_plot() +
    ggplot2::theme(
      panel.border = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank()
    )
}
