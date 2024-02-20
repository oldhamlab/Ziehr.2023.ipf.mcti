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

read_mims <- function(file) {
  img <- mims::mims_read(file)

  n_ratio <- mims::mims_ratio(img, "12C15N", "12C14N")
  h_ratio <- mims::mims_ratio(img, "12C22H", "12C21H")

  idx <- which(img@metadata$globalMetadata$symbols == "12C14N")
  mims::mims_sum_frames(img, idx) |>
    EBImage::abind(n_ratio[,,1], along = 3) |>
    EBImage::abind(h_ratio[,,1], along = 3)
}

extract_ratios <- function(img, nm) {
  histo <- EBImage::normalize(img[,,1])

  threshold <- EBImage::otsu(histo)
  mask <- histo > threshold

  # EBImage::display(mask)

  n_ratio <- img[,,2][mask]
  h_ratio <- img[,,3][mask]
  tibble::tibble(
    name = nm,
    n_ratio = n_ratio,
    h_ratio = h_ratio
  )
}

format_ratios <- function(x) {
  dplyr::bind_rows(x) |>
    tidyr::separate(name, c("group", "region", "replicate"), sep = "_") |>
    dplyr::mutate(
      n_ratio_norm = n_ratio / 0.0037,
      h_ratio_norm = h_ratio / 0.0001,
      group = toupper(group),
      region = factor(
        region,
        levels = c("alv", "fib"),
        labels = c("Alveolar", "Fibrotic")
      )
    ) |>
    refactor() |>
    dplyr::mutate(group = forcats::fct_drop(group))
}
