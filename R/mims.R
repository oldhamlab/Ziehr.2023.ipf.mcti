# mims.R

clean_mims <- function(files) {
  names(files) <- stringr::str_replace(basename(files), "(-interleaved)*_data.csv", "")

  df_slice <-
    purrr::map(
      files[!stringr::str_detect(files, "interleaved")],
      \(x) readr::read_csv(x, show_col_types = FALSE)
    ) |>
    purrr::list_rbind(names_to = "file")

  purrr::map(
    files[stringr::str_detect(files, "interleaved")],
    \(x) readr::read_csv(x, show_col_types = FALSE)
  ) |>
    purrr::list_rbind(names_to = "file") |>
    dplyr::left_join(
      df_slice,
      by = c("file", "Roi group", "Roi tags", "Roi name")
    ) |>
    dplyr::rename(
      feature = "Roi group",
      tag = "Roi tags",
      roi_id = "Roi name"
    ) |>
    tidyr::pivot_longer(
      cols = tidyselect::matches("\\d+\\w+"),
      names_to = "measurement",
      values_to = "value"
    ) |>
    tidyr::separate_wider_delim("measurement", names = c("ion", "measurement"), delim = " ") |>
    tidyr::pivot_wider(
      names_from = "measurement",
      values_from = "value"
    ) |>
    tidyr::separate_wider_delim(
      cols = file,
      names = c("group", "region", "image", "date", "sample"),
      delim = "_"
    ) |>
    dplyr::mutate(group = toupper(group)) |>
    refactor()
}
