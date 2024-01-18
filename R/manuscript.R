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
