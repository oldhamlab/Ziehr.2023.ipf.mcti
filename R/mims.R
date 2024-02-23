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
  symbols <- img@metadata$globalMetadata$symbols
  empty <- EBImage::Image(NA, dim(img)[1:2])
  n_ratio <- h_ratio <- c_ratio <- empty

  if ("12C15N" %in% symbols) {
    n_ratio <-
      mims::mims_ratio(img, "12C15N", "12C14N")[,,1]
  }
  if ("12C22H" %in% symbols) {
    h_ratio <-
      mims::mims_ratio(img, "12C22H", "12C21H")[,,1]
  }
  if ("12C13C" %in% symbols) {
    c_ratio <-
      mims::mims_ratio(img, "12C13C", "12C2")[,,1]
  }

  idx <- which(symbols == "12C14N")
  mims::mims_sum_frames(img, idx) |>
    EBImage::abind(n_ratio, along = 3) |>
    EBImage::abind(h_ratio, along = 3) |>
    EBImage::abind(c_ratio, along = 3)
}

extract_ratios <- function(img, nm) {
  histo <- EBImage::normalize(img[,,1])

  threshold <- EBImage::otsu(histo)
  mask <- histo > threshold

  # EBImage::display(mask)

  n_ratio <- img[,,2][mask]
  h_ratio <- img[,,3][mask]
  c_ratio <- img[,,4][mask]
  tibble::tibble(
    name = nm,
    n_ratio = n_ratio,
    h_ratio = h_ratio,
    c_ratio = c_ratio
  )
}

format_ratios <- function(x) {
  x |>
    tidyr::separate(
      name,
      c("group", "animal", "region", "replicate"),
      sep = "_"
    ) |>
    dplyr::mutate(
      n_ratio_norm = n_ratio / (0.368/99.632),
      h_ratio_norm = h_ratio / (0.0115/99.9885),
      c_ratio_norm = c_ratio / (1.07/98.93),
      group = toupper(group),
      region = factor(
        region,
        levels = c("alv", "fib"),
        labels = c("Alveolar", "Fibrotic")
      )
    ) |>
    tidyr::pivot_longer(
      -c("group":"replicate"),
      names_to = c("element", ".value"),
      names_pattern = c("(\\w)_(.+)"),
      cols_vary = "slowest"
    ) |>
    dplyr::filter(!is.na(.data$ratio)) |>
    dplyr::filter(.data$ratio > 0) |>
    dplyr::mutate(
      tracer = dplyr::case_when(
        element == "n" ~ "proline",
        element == "h" & animal != 1 ~ "lactate",
        element == "h" & animal == 1 ~ "glucose",
        element == "c" ~ "glucose"
      ),
      .after = "element"
    ) |>
    refactor() |>
    dplyr::mutate(group = forcats::fct_drop(group))
}

average_mims <- function(x) {
  x |>
    dplyr::group_by(dplyr::across(c("group":"tracer"))) |>
    dplyr::summarise(
      count = dplyr::n(),
      ratio_mean = mean(ratio, trim = 0.1),
      ratio_sd = sd(ratio),
      ratio_norm_mean = mean(ratio_norm, trim = 0.1),
      ratio_norm_sd = sd(ratio_norm)
    ) |>
    refactor()
}

summarize_mims <- function(x) {
  df <-
    x |>
    dplyr::filter(tracer != "lactate") |>
    dplyr::group_by(group, animal, tracer) |>
    dplyr::summarise(value = sum(count * ratio_norm_mean) / sum(count)) |>
    dplyr::group_by(tracer) |>
    dplyr::mutate(
      tracer = factor(
        tracer,
        levels = c("proline", "glucose"),
        labels = c("Proline", "Glucose")
      ),
      measurement = "tracing",

      value = (value - 1) / mean(value[group == "Ctl"] - 1)
    ) |>
    normalize("value", "animal") |>
    dplyr::arrange(tracer, group, animal)
}

stats_mims <- function(x) {
  x |>
    dplyr::group_by(tracer) |>
    tidyr::nest() |>
    dplyr::mutate(stats = purrr::map(
      data,
      \(x) stats_histo(x, measure = "tracing", mixed = TRUE)
    )) |>
    tidyr::unnest(c(stats)) |>
    refactor()
}
