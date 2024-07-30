# blots.R

clean_blots <- function(df) {
  new_names <-
    c(
      "sma" = "α-SMA",
      "pampk" = "pAMPK",
      "ampk" = "AMPK",
      "psmad3" = "pSmad3",
      "smad3" = "Smad3",
      "erk" = "ERK",
      "perk" = "pERK",
      "lactylk" = "Kla",
      "mct1" = "MCT1",
      "mct4" = "MCT4",
      "col1a1" = "Col1a1",
      "col1a2" = "Col1a2",
      "fn1" = "FN1",
      "col3a1" = "Col3a1",
      "cnn1" = "CNN1",
      "hif1a" = "HIF-1α",
      "akt" = "AKT",
      "pakt" = "pAKT",
      "h3" = "H3"
    )

  df |>
    dplyr::group_by(dplyr::across(c("experiment":"trial", "protein"))) |>
    dplyr::mutate(
      dplyr::across(
        "background":"density",
        \(x) x / mean(x, na.rm = TRUE)
      ),
      ratio = .data$density / .data$background
    ) |>
    dplyr::group_by(dplyr::across("experiment":"trial")) |>
    dplyr::mutate(protein = new_names[.data$protein]) |>
    refactor() |>
    dplyr::mutate(group = dplyr::cur_group_id()) |>
    dplyr::ungroup()
}

filter_blots <- function(df) {
  prots <- c("AKT", "pAKT", "AMPK", "pAMPK", "ERK", "pERK")
  df |>
    dplyr::filter(!(.data$experiment == "ipf" & .data$batch == 2)) |>
    dplyr::filter(!(.data$experiment == "ipf-lf" & .data$batch %in% 1:2)) |>
    dplyr::filter(!(.data$experiment == "dual" & .data$protein == "Col1a1" & .data$batch %in% 1:4)) |>
    # dplyr::filter(!(.data$experiment == "dual" & .data$protein %in% prots & .data$batch %in% 1:4)) |>
    dplyr::filter(!(.data$experiment == "dual" & .data$protein %in% c("pAKT", "AKT") & .data$batch == 5 & .data$trial == 2)) |>
    dplyr::filter(!(.data$experiment == "dual" & .data$protein %in% c("pERK", "ERK") & .data$batch == 3 & .data$trial %in% 3:4)) |>
    dplyr::filter(!(.data$experiment == "six" & .data$protein == "HIF-1α" & .data$trial == 5))
}

index_blots <- function(df) {
  df |>
    dplyr::select("experiment", "protein") |>
    dplyr::distinct() |>
    dplyr::mutate(
      names = stringr::str_c(.data$experiment, .data$protein, sep = "_")
    ) |>
    dplyr::arrange(.data$names)
}

norm_blot <- function(df, exp, prot, treat = NULL, cond = NULL) {
  if (is.null(treat)) {
    treat <- unique(df$treatment)
  }
  if (is.null(cond)) {
    cond <- unique(df$condition)
  }

  df |>
    dplyr::group_by(.data$protein) |>
    dplyr::filter(
      .data$experiment %in% exp &
        .data$protein %in% prot &
        .data$treatment %in% treat &
        .data$condition %in% cond
    ) |>
    dplyr::mutate(
      norm = .data$ratio / mean(
        .data$ratio[
          .data$condition == min(.data$condition) &
            .data$treatment == min(.data$treatment)
        ],
        na.rm = TRUE
      ),
      protein = factor(.data$protein, levels = prot)
    )
}

analyze_blot <- function(df, paired, comp = NULL) {
  # require multiple groups
  if (length(unique(df$group)) == 1 &&
      all(df$experiment %nin% c("comp", "ipf"))) {
    return()
  }

  if (length(unique(df$treatment)) == 1) {
    out <- ratio_ttest(df, "norm")
  } else {
    out <- twofactor(df, "norm", mixed = paired, comp)
  }
  out |>
    dplyr::mutate(protein = unique(df$protein), .before = 1)
}

plot_blot <- function(df, stats = NULL, title = NULL) {
  if (length(unique(df$treatment)) == 1) {
    plot_one_factor(df, stats, "condition", "norm", title, "Fold change")
  } else {
    plot_two_factor(df, stats, "treatment", "norm", title, "Fold change")
  }
}

blot_ratios <- function(df) {
  num <- c("pAKT", "pAMPK", "pERK", "pSmad3", "Kla")
  denom <- c("AKT", "AMPK", "ERK", "Smad3", "H3")
  ratios <- stringr::str_c(num, denom, sep = "/")

  df |>
    dplyr::filter(
      stringr::str_detect(.data$protein, paste(c(num, denom), collapse = "|"))
    ) |>
    dplyr::mutate(
      temp = ifelse(
        .data$protein %in% num,
        "density",
        "background"
      ),
      protein = purrr::map_chr(
        .data$protein,
        \(x) ratios[stringr::str_detect(ratios, x)]
      ),
      .after = "dose"
    ) |>
    dplyr::select("group", "experiment":"protein", "ratio") |>
    tidyr::pivot_wider(
      names_from = "temp",
      values_from = "ratio"
    ) |>
    dplyr::mutate(ratio = .data$density / .data$background) |>
    dplyr::filter(!is.na(.data$ratio))
}
