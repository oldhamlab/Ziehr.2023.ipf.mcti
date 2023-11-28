# utils.R

# I/O ---------------------------------------------------------------------

list_files <- function(nm, path) {
  list.files(
    path = path,
    pattern = nm,
    all.files = TRUE,
    full.names = TRUE,
    recursive = TRUE,
    include.dirs = TRUE
  )
}

raw_data_path <- function(nm) {
  list_files(nm = nm, path = "data-raw")
}

manuscript_path <- function(nm) {
  list_files(nm = nm, path = "manuscript")
}

plot_path <- function(nm) {
  path <- stringr::str_c("analysis/figures/", nm)
  if (dir.exists(path)) unlink(path, recursive = TRUE)
  if (!dir.exists(path)) dir.create(path = path, recursive = TRUE)
  path
}


# helpers -----------------------------------------------------------------

"%nin%" <- function(x, table) {
  match(x, table, nomatch = 0L) == 0L
}


# factors -----------------------------------------------------------------

refactor <- function(df) {
  if ("condition" %in% names(df)) {
    cond <- c(
      "Ctl" = "control",
      "Ctl" = "Ctl",
      "TGFβ" = "tgfb",
      "TGFβ" = "TGFb",
      "TGFβ" = "TGFβ",
      "IPF" = "ipf",
      "Bleo" = "bleo",
      "Ctl-LF" = "nhlf",
      "IPF-LF" = "ipf-lf",
      NULL
    )
    df <-
      dplyr::mutate(
        df,
        condition = factor(
          .data$condition,
          levels = cond,
          labels = names(cond),
          ordered = TRUE
        )
      )
  }
  if ("treatment" %in% names(df)) {
    treat <- c(
      "None" = "none",
      "Veh" = "DMSO",
      "Veh" = "vehicle",
      "Veh" = "Veh",
      "AZD" = "AZD3965",
      "AZD" = "AZD",
      "AR" = "AR-C155858",
      "VB" = "VB124",
      "VB" = "VB",
      "AZD/VB" = "Dual",
      "AZD/VB" =  "AZD/VB",
      "AR/VB" = "AR/VB",
      "DIC" = "diclofenac",
      "KET" = "ketoprofen",
      "SA" = "salicylate",
      "siCTL" = "siCTL",
      "siMCT1" = "siMCT1",
      "siMCT4" = "siMCT4",
      "siMCT1/4" = "siMCT1/4",
      "Lac" = "lactate",
      NULL
    )
    df <-
      dplyr::mutate(
        df,
        treatment = factor(
          .data$treatment,
          levels = treat,
          labels = names(treat),
          ordered = TRUE
        )
      )
  }

  if ("group" %in% names(df)) {
    grp <- c(
      "Ctl\nNone" = "Ctl\nNone",
      "TGFβ\nNone" = "TGFβ\nNone",
      "Ctl\nVeh" = "Ctl\nVeh",
      "TGFβ\nVeh" = "TGFβ\nVeh",
      "Ctl\nAZD" = "Ctl\nAZD",
      "TGFβ\nAZD" = "TGFβ\nAZD",
      "Ctl\nVB" = "Ctl\nVB",
      "TGFβ\nVB" = "TGFβ\nVB",
      "Ctl\nAZD/VB" = "Ctl\nAZD/VB",
      "TGFβ\nAZD/VB" = "TGFβ\nAZD/VB",
      "Ctl" = "Ctl",
      "Ctl" = "CONTROL",
      "TGFβ" = "TGFβ",
      "IPF" = "IPF",
      "Bleo" = "Bleo",
      "Bleo" = "BLEO",
      "Ctl\nVeh" = "Ctl_Veh",
      "TGFβ\nVeh" = "TGFβ_Veh",
      "TGFβ\nAZD" = "TGFβ_AZD",
      "TGFβ\nVB" = "TGFβ_VB",
      "TGFβ\nAZD/VB" = "TGFβ_AZD/VB",
      "Ctl\nNone" = "control\nnone",
      "Ctl\nVeh" = "control\nDMSO",
      "Ctl\nAZD" = "control\nAZD3965",
      "AZD" = "AZD",
      "Ctl\nVB" = "control\nVB124",
      "VB" = "VB",
      "Ctl\nAZD/VB" = "control\nAZD/VB",
      "TGFβ\nNone" = "tgfb\nnone",
      "TGFβ\nVeh" = "tgfb\nDMSO",
      "TGFβ\nAZD" = "tgfb\nAZD3965",
      "TGFβ\nVB" = "tgfb\nVB124",
      "TGFβ\nAZD/VB" = "tgfb\nAZD/VB" ,
      NULL
    )
    df <-
      dplyr::mutate(
        df,
        group = factor(
          .data$group,
          levels = grp,
          labels = names(grp)
        )
      )
  }
  df
}


# interpolation -----------------------------------------------------------

make_std_curves <- function(df, fo = NULL) {
  if (is.null(fo)){
    fo <- \(x) stats::lm(value ~ conc, data = x, na.action = modelr::na.warn)
  }

  df |>
    dplyr::filter(!is.na(.data$conc)) |>
    dplyr::select(!tidyselect::where(\(x) all(is.na(x)))) |>
    dplyr::group_by(dplyr::across(-c("conc", "value"))) |>
    tidyr::nest() |>
    {\(x) dplyr::mutate(
      x,
      title = stringr::str_c(!!!rlang::syms(dplyr::groups(x)), sep = "_")
    )}() |>
    dplyr::ungroup() |>
    dplyr::mutate(
      model = furrr::future_map(.data$data, fo),
      summary = furrr::future_map(.data$model, \(x) broom::glance(x)),
      plots = furrr::future_map2(.data$data, .data$title, make_std_plots)
    ) |>
    dplyr::group_by(
      dplyr::across(
        -c("data", "title", "model", "summary", "plots")
      )
    )
}

make_std_plots <- function(df, title = NULL) {
  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = .data$conc,
      y = .data$value
    ) +
    ggplot2::geom_smooth(
      method = stats::lm,
      formula = y ~ x,
      color = "gray20",
      se = FALSE
    ) +
    ggplot2::geom_point(
      size = 3,
      color = "darkblue"
    ) +
    ggplot2::stat_summary(
      fun = "mean",
      size = 4,
      geom = "point",
      alpha = 0.8,
      color = "darkblue"
    ) +
    ggplot2::labs(
      x = "Concentration",
      y = "Value",
      title = title
    )
}

write_plot <- function(p, nm, path, width = 15, height = 20, units = "in", ...) {
  nm <- stringr::str_c(nm, ".png")
  ggplot2::ggsave(
    plot = p,
    filename = nm,
    path = path,
    device = ragg::agg_png,
    res = 300,
    ...
  )
  invisible(stringr::str_c(path, "/", nm))
}

write_plot_list <- function(p_list, nm_list, path, ...) {
  path <- plot_path(path)
  furrr::future_walk2(
    p_list,
    nm_list,
    \(x, y) write_plot(x, y, path)
  )
  invisible(path)
}

interp_data <- function(df, std) {
  df |>
    dplyr::filter(is.na(.data$conc)) |>
    dplyr::select(-"conc") |>
    dplyr::group_by(dplyr::across(dplyr::group_vars(std))) |>
    tidyr::nest() |>

    ### PROBLEMS WITH LEFT JOIN HERE ###
    dplyr::left_join(dplyr::select(std, "model")) |>
    dplyr::mutate(
      conc = purrr::map2(
        .data$data,
        .data$model, wmo::interpolate
      )
    ) |>
    tidyr::unnest(c("data", "conc")) |>
    dplyr::select(-c("model", "value"))
}
