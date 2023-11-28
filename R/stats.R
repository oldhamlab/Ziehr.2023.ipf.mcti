# stats.R

annot_p <- function(vec) {
  dplyr::case_when(
    vec < 0.05 ~ "*",
    .default = ""
  )
}


annot_stats <- function(stats) {
  col <- intersect(names(stats), c("adj.p.value", "p.value"))
  stats |>
    dplyr::mutate(
      y = Inf,
      vjust = 1.5,
      label = annot_p(.data[[col]])
    )
}


ratio_ttest <- function(df, y, paired){
  fo <- stats::formula(paste("log(", y, ") ~ condition"))

  stats::t.test(fo, data = df, paired = paired) |>
    broom::tidy() |>
    dplyr::select("p.value") |>
    dplyr::mutate(
      condition = NA_character_,
      x = 1.5
    ) |>
    annot_stats()
}


onefactor <- function(
    df,
    y,
    mixed = TRUE,
    comps,
    fo = NULL,
    emm = NULL
) {

  x <-
    df |>
    dplyr::mutate(
      dplyr::across(c("condition", "treatment"), forcats::fct_drop)
    )

  withCallingHandlers(
    message = \(cnd) rlang::warn(cnd$message),
    {
      if (mixed) {
        fo <- stats::formula(paste(y, " ~ treatment + (1 | group)"))
        m <- lmerTest::lmer(fo, data = x)
      } else {
        fo <- stats::formula(paste(y, " ~ treatment"))
        m <- stats::lm(fo, data = x)
      }
    }
  )

  if (is.null(emm)) {
    emm = stats::formula("~ treatment")
  }

  m |>
    emmeans::emmeans(emm) |>
    emmeans::contrast(method = comps, adjust = "dunnettx") |>
    broom::tidy() |>
    dplyr::select("contrast", tidyselect::contains("p.value")) |>
    tidyr::separate(
      .data$contrast,
      c("treatment", "condition"),
      sep = "\\.",
      fill = "right"
    ) |>
    dplyr::mutate(
      treatment = factor(.data$treatment, levels = levels(x$treatment)),
      condition = factor(.data$condition, levels = levels(x$condition))
    ) |>
    tidyr::complete(.data$treatment) |>
    annot_stats()
}


twofactor <- function(df, y, mixed = TRUE, comps, log = FALSE) {
  x <-
    df |>
    dplyr::mutate(
      dplyr::across(c("condition", "treatment"), forcats::fct_drop)
    )

  if (log) {
    y <- paste0("log(", y, ")")
  }

  withCallingHandlers(
    message = \(cnd) rlang::warn(cnd$message),
    {
      if (mixed) {
        fo <- stats::formula(paste(y, " ~ condition * treatment + (1 | group)"))
        m <- lmerTest::lmer(fo, data = x)
      } else {
        fo <- stats::formula(paste(y, " ~ condition * treatment"))
        m <- stats::lm(fo, data = x)
      }
    }
  )

  m |>
    emmeans::emmeans(~ condition * treatment) |>
    emmeans::contrast(method = comps, adjust = "mvt") |>
    broom::tidy() |>
    dplyr::select("contrast", tidyselect::contains("p.value")) |>
    tidyr::separate(
      .data$contrast,
      c("treatment", "condition"),
      sep = "\\.",
      fill = "right"
    ) |>
    dplyr::mutate(
      treatment = factor(.data$treatment, levels = levels(x$treatment)),
      condition = factor(.data$condition, levels = levels(x$condition))
    ) |>
    tidyr::complete(.data$condition, .data$treatment) |>
    annot_stats()
}

group_twofactor <- function(df, y, comps, ...) {
  cols <- c("protein", "name", "rate", "stage", "stats", "measurement")
  df |>
    tidyr::nest() |>
    dplyr::mutate(stats = purrr::map(
      .data$data,
      twofactor,
      y = y,
      comps = comps,
      ...
    )) |>
    dplyr::select(tidyselect::any_of(cols)) |>
    tidyr::unnest(c("stats"))
}
