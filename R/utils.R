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
