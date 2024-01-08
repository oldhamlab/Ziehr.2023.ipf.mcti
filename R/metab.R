# metab.R


# utils -------------------------------------------------------------------

read_pathways <- function(filename){
  cols <- c("entrez_gene_ids", "metabolites")
  readr::read_tsv(filename, col_types = "cccc") |>
    tidyr::unite("pathway", source, pathway, sep = " | ") |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(cols),
        \(x) stringr::str_split(x, pattern = ",")
      )
    ) |>
    dplyr::select(pathway, tidyselect::any_of(cols)) |>
    tibble::deframe()
}

new_tbl_se <- function(
    tbl,
    a_data,
    f_names,
    f_data = NULL,
    s_names,
    s_data = NULL
){
  structure(
    tbl,
    class = c("tbl_se", class(tbl)),
    a_data = a_data,
    f_names = f_names,
    s_names = s_names,
    f_data = f_data,
    s_data = s_data
  )
}

tbl_to_se <- function(tbl_se, assay_name){
  assay_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "s_names"),
      attr(tbl_se, "a_data")
    ) |>
    tidyr::pivot_wider(
      names_from = attr(tbl_se, "s_names"),
      values_from = attr(tbl_se, "a_data")
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names"))

  feature_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "f_names"),
      attr(tbl_se, "f_data")
    ) |>
    dplyr::group_by(!!rlang::sym(attr(tbl_se, "f_names"))) |>
    dplyr::summarise(
      metabolite = unique(metabolite),
      mz = mean(mz, na.rm = TRUE),
      mz_min = min(mz, na.rm = TRUE),
      mz_max = max(mz, na.rm = TRUE),
      rt = mean(rt, na.rm = TRUE),
      rt_min = min(rt, na.rm = TRUE),
      rt_max = max(rt, na.rm = TRUE)
    ) |>
    tibble::column_to_rownames(attr(tbl_se, "f_names")) |>
    {\(x) x[match(rownames(assay_data), rownames(x)), ]}()

  sample_data <-
    tbl_se |>
    dplyr::select(
      attr(tbl_se, "s_names"),
      attr(tbl_se, "s_data")
    ) |>
    dplyr::distinct() |>
    tibble::column_to_rownames(attr(tbl_se, "s_names")) |>
    {\(x) x[match(colnames(assay_data), rownames(x)), ]}()

  SummarizedExperiment::SummarizedExperiment(
    assays = assay_data,
    rowData = feature_data,
    colData = sample_data
  )
}

se_to_tibble <- function(se){
  tibble::as_tibble(SummarizedExperiment::assay(se), rownames = "feature") |>
    tidyr::pivot_longer(-feature, names_to = "sample", values_to = "area") |>
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::colData(se), rownames = "sample"),
      by = "sample"
    ) |>
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(se), rownames = "feature"),
      by = "feature"
    )
}

remove_missing_metab <- function(raw) {
  qc <- SummarizedExperiment::assay(raw[, raw$type == "qc"])
  missing <- names(which(apply(qc, 1, function(x) sum(is.na(x))) > 0))
  raw[rownames(raw) %nin% missing, raw$type %nin% c("water", "blank")]
}

correct_drift <- function(missing) {
  prepare_assay_data <- function(se)  {
    SummarizedExperiment::assay(se) |>
      tibble::rownames_to_column("hmdb") |>
      tidyr::pivot_longer(-"hmdb", names_to = "sample", values_to = "value") |>
      dplyr::mutate(
        value = log(.data$value),
        run_order = as.numeric(stringr::str_extract(.data$sample, "\\d{2}"))
      ) |>
      dplyr::group_by(.data$hmdb) |>
      tidyr::nest()
  }

  models <-
    missing[, missing$type == "qc"] |>
    prepare_assay_data() |>
    dplyr::mutate(
      model = map(data, \(x) smooth.spline(x = x$run_order, y = x$value, spar = 0.2)),
      mean = map_dbl(data, \(x) mean(x$value))
    ) |>
    dplyr::select(-"data")

  corrected <-
    missing |>
    prepare_assay_data() |>
    dplyr::left_join(models, by = "hmdb") |>
    dplyr::mutate(pred = purrr::map2(model, data, \(x, y) predict(x, y$run_order)$y)) |>
    tidyr::unnest(c("data", "pred")) |>
    dplyr::mutate(corr = value + mean - pred) |>
    dplyr::select("hmdb", "sample", "corr") |>
    tidyr::pivot_wider(names_from = "sample", values_from = "corr") |>
    tibble::column_to_rownames("hmdb") |>
    exp()

  SummarizedExperiment::assay(missing) <- corrected
  missing
}

quality_control <- function(drift) {
  rsd <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, \(x) 1.4826 * mad(x, na.rm = TRUE) / median(x, na.rm = TRUE))

  mad_qc <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, \(x) mad(x, na.rm = TRUE))

  ref <-
    drift[, drift$type == "qc"] |>
    SummarizedExperiment::assay() |>
    apply(1, \(x) median(x, na.rm = TRUE))

  mad_s <-
    drift[, drift$type == "sample"] |>
    SummarizedExperiment::assay() |>
    apply(1, \(x) mad(x, na.rm = TRUE))

  d_ratio <- mad_qc / mad_s

  SummarizedExperiment::rowData(drift)$rsd <- rsd
  SummarizedExperiment::rowData(drift)$rsd <- d_ratio
  SummarizedExperiment::rowData(drift)$good <- rsd < 0.2 & d_ratio < 0.4
  SummarizedExperiment::rowData(drift)$reference <- ref

  drift[SummarizedExperiment::rowData(drift)$good == TRUE, drift$type != "qc"]
}

impute_missing <- function(qc) {
  set.seed(42)
  SummarizedExperiment::assay(qc) <-
    missForest::missForest(
      t(SummarizedExperiment::assay(qc)),
      maxiter = 10
    )$ximp |>
    t()
  qc
}

pqn <- function(imputed) {
  mat <- SummarizedExperiment::assay(imputed)
  quotients <- mat / SummarizedExperiment::rowData(imputed)$reference
  quotient_medians <- apply(quotients, 2, median)
  SummarizedExperiment::assay(imputed) <- t(t(mat) / quotient_medians)
  imputed
}

annot_metabs <- function(se) {
  fdata <-
    SummarizedExperiment::rowData(se) |>
    data.frame() |>
    tibble::as_tibble(rownames = "HMDB")

  ah <- AnnotationHub::AnnotationHub()
  df <-
    ah[["AH91792"]] |>
    dplyr::filter(HMDB %in% fdata$HMDB) |>
    dplyr::group_by(HMDB) |>
    dplyr::arrange(HMDB, KEGG, ChEBI, .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::select(-"Name")

  annot_fdata <-
    dplyr::left_join(fdata, df, by ="HMDB") |>
    tibble::column_to_rownames("HMDB")

  SummarizedExperiment::rowData(se) <- annot_fdata
  se
}

plot_pca_metab <- function(clean, batch = TRUE, show_reps = TRUE){
  df <-
    SummarizedExperiment::assay(clean) |>
    log() |>
    t() |>
    scale() |>
    t()

  pheno <-
    SummarizedExperiment::colData(clean) |>
    tibble::as_tibble() |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.factor), forcats::fct_drop))

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- stringr::str_replace(colnames(design), "group", "")

  if (batch) {
    df <-
      limma::removeBatchEffect(df, batch = pheno$replicate, design = design)
  }

  df <-
    df |>
    t() |>
    pcaMethods::pca(scale = "none", center = TRUE)

  percent_variance <- round(100 * c(df@R2[[1]], df@R2[[2]]))

  p <-
    merge(
      pcaMethods::scores(df),
      SummarizedExperiment::colData(clean),
      by = 0
    ) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = group
    ) +
    ggforce::geom_mark_ellipse(
      linewidth = 0.25,
      expand = ggplot2::unit(2, "mm"),
      label.fontsize = 5,
      label.fontface = "plain",
      label.family = "Calibri",
      label.hjust = 0.5,
      label.buffer = ggplot2::unit(0, "mm"),
      label.margin = ggplot2::margin(-1.5, -1.5, -1.5, -1.5, "mm"),
      con.type = "none",
      show.legend = FALSE
    ) +
    ggplot2::labs(
      x = paste0("PC1: ", percent_variance[1], "% variance"),
      y = paste0("PC2: ", percent_variance[2], "% variance"),
      shape = NULL
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      breaks = c("TGFβ\nNone", "TGFβ\nVeh", "TGFβ\nAZD", "TGFβ\nVB", "TGFβ\nAZD/VB", "Ctl", "Bleo", "AZD", "VB"),
      labels = c("None", "Veh", "AZD", "VB", "AZD/VB", "Ctl", "Bleo", "AZD", "VB"),
      values = clrs,
      aesthetics = c("color", "fill"),
      limits = force
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0.1, 0.05)),
      limits = nice_limits
    ) +
    ggplot2::guides(shape = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::coord_fixed(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.height = ggplot2::unit(1, "lines"),
      legend.box = "vertical"
    ) +
    NULL

  if (show_reps) {
    p +
      ggplot2::geom_point(
        ggplot2::aes(
          fill = group,
          shape = replicate
        ),
        size = 2,
        show.legend = TRUE
      )
  } else {
    p +
      ggplot2::geom_point(
        ggplot2::aes(
          fill = group,
        ),
        pch = 21,
        color = "white",
        size = 1.5,
        show.legend = TRUE
      )
  }
}

fit_limma_metab <- function(se, batch = TRUE, ...){
  input <-
    SummarizedExperiment::assay(se) |>
    log(base = 2)

  pheno <-
    SummarizedExperiment::colData(se, ...) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      group = forcats::fct_relabel(group, stringr::str_replace_all, "(\n|/)", "."),
      group = forcats::fct_drop(group)
    )

  design <- model.matrix(~ 0 + group, data = pheno)
  colnames(design) <- stringr::str_replace(colnames(design), "group", "")

  cm <-
    limma::makeContrasts(
      ...,
      levels = design
    )

  if (batch) {
    corfit <- limma::duplicateCorrelation(input, design, block = pheno$replicate)
    fit <-
      limma::lmFit(
        input,
        design,
        block = pheno$replicate,
        correlation = corfit$consensus
      )
  } else {
    fit <-
      limma::lmFit(
        input,
        design
      )
  }

  fit |>
    limma::contrasts.fit(cm) |>
    limma::eBayes()
}

top_table_metab <- function(se, fit, contrast){
  limma::topTable(
    fit,
    number = Inf,
    p.value = 1,
    coef = contrast
  ) |>
    tibble::as_tibble(rownames = "hmdb") |>
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(se), rownames = "hmdb"),
      by = "hmdb"
    )
}

plot_mois <- function(se, tt, moi){
  df <-
    se_to_tibble(se) |>
    dplyr::filter(metabolite %in% moi) |>
    dplyr::group_by(metabolite) |>
    dplyr::mutate(
      area = area / mean(area[condition == min(condition) & treatment == min(treatment)]),
      metabolite = dplyr::case_when(
        metabolite == "2-hydroxyglutarate" ~ "2HG",
        metabolite == "glyceraldehyde 3-phosphate" ~ "GAP",
        TRUE ~ metabolite
      ),
      metabolite = factor(metabolite, levels = moi)
    )

  stats <-
    tt |>
    dplyr::filter(metabolite %in% moi) |>
    dplyr::select(treatment, metabolite, adj.P.Val) |>
    dplyr::mutate(
      y = Inf,
      vjust = 2,
      label = annot_p(adj.P.Val),
      treatment = factor(treatment, levels = c("Veh", "AZD", "VB", "AZD/VB"))
    )

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(
      ggplot2::vars(.data$metabolite),
      labeller = ggplot2::as_labeller(toupper),
      nrow = 1
    ) +
    ggplot2::aes(
      x = .data$treatment,
      y = .data$area,
      fill = .data$group
    ) +
    ggplot2::geom_bar(
      stat = "summary",
      fun = "mean",
      ggplot2::aes(fill = .data[["group"]]),
      position = ggplot2::position_dodge(width = 0.55),
      width = 0.5,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(group = .data[["condition"]]),
      position = ggplot2::position_dodge(width = 0.55),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.25,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      ggplot2::aes(group = .data[["condition"]]),
      dodge.width = 0.55,
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = stats,
      ggplot2::aes(
        x = .data$treatment,
        y = .data$y,
        label = .data$label,
        vjust = .data$vjust
      ),
      color = "black",
      inherit.aes = FALSE
    ) +
    ggplot2::scale_fill_manual(values = clrs) +
    ggplot2::labs(
      x = NULL,
      y = "Peak area (normalized)",
      fill = NULL
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}

plot_vol_metab <- function(
    tt,
    fills = c("Ctl", "IPF"),
    title = NULL,
    mois = NULL
) {
  left <-
    tt |>
    dplyr::filter(logFC < 0 & adj.P.Val < 0.1) |>
    dplyr::slice_min(t, n = 10)

  right <-
    tt |>
    dplyr::filter(logFC > 0 & adj.P.Val < 0.1) |>
    dplyr::slice_max(t, n = 10)

  nudge <- max(0.8 * abs(c(min(tt$logFC), max(tt$logFC))))
  left_nudge <- floor(min(tt$logFC) - nudge)
  right_nudge <- ceiling(max(tt$logFC) + nudge)

  ggplot2::ggplot(tt) +
    ggplot2::aes(
      x = .data$logFC,
      y = .data$adj.P.Val
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% mois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      xlim = c(-Inf, Inf),
      nudge_x = left_nudge - left$logFC,
      hjust = 0,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = metabolite,
        color = metabolite %in% mois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      xlim = c(-Inf, Inf),
      nudge_x = right_nudge - right$logFC,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(tt, adj.P.Val > 0.1),
      pch = 21,
      color = "white",
      fill = "grey80"
    ) +
    ggplot2::geom_point(
      data = subset(tt, logFC > 0 & adj.P.Val < 0.1),
      pch = 21,
      color = "white",
      fill = clrs[[fills[[2]]]]
    ) +
    ggplot2::geom_point(
      data = subset(tt, logFC < 0 & adj.P.Val < 0.1),
      pch = 21,
      color = "white",
      fill = clrs[[fills[[1]]]]
    ) +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log(),
      breaks = scales::breaks_log(6)
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 7),
      # expand = ggplot2::expansion(mult = 0.4),
      labels = scales::label_math(2^.x)
    ) +
    ggplot2::labs(
      x = "Fold change",
      y = "Adjusted P value",
      title = title
    ) +
    theme_plot() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.y = ggplot2::element_text(hjust = 0)
    ) +
    NULL
}

read_pathways <- function(filename){
  cols <- c("entrez_gene_ids", "metabolites")
  readr::read_tsv(filename, col_types = "cccc") |>
    tidyr::unite("pathway", source, pathway, sep = " | ") |>
    dplyr::mutate(
      dplyr::across(
        tidyselect::any_of(cols),
        \(x) stringr::str_split(x, pattern = ",")
      )
    ) |>
    dplyr::select(pathway, tidyselect::any_of(cols)) |>
    tibble::deframe()
}

run_msea <- function(tt, pathways){
  stats <-
    tt |>
    dplyr::filter(!is.na(KEGG)) |>
    dplyr::mutate(KEGG = stringr::str_c("kegg:", KEGG)) |>
    dplyr::select(KEGG, t) |>
    tibble::deframe()

  fgsea::fgsea(
    pathways = pathways,
    stats = stats,
    minSize = 3,
    BPPARAM = BiocParallel::bpparam()
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(pval < 0.1) |>
    tidyr::separate(pathway, c("source", "pathway"), sep = " \\| ") |>
    dplyr::arrange(desc(NES))
}

plot_msea <- function(df, src = "KEGG", title = NULL, colors = NULL){
  x <-
    df |>
    dplyr::filter(source %in% src) |>
    dplyr::mutate(
      pathway = stringr::str_replace(pathway, " - Homo sapiens \\(human\\)", ""),
      pathway = forcats::fct_reorder(pathway, NES)
    ) |>
    dplyr::arrange(NES)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = NES,
      y = pathway,
      fill = NES >= 0
    ) +
    ggplot2::geom_col(
      show.legend = FALSE
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.25
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits
    ) +
    ggplot2::scale_y_discrete(position = "right") +
    ggplot2::scale_fill_manual(values = unname(clrs[colors])) +
    ggplot2::labs(
      y = NULL,
      title = title,
      x = "Normalized enrichment score"
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5)
    ) +
    NULL
}

plot_msea_table <- function(df, title, clr) {
  clr <- clrs[clr]
  x <-
    df |>
    dplyr::filter(.data$source == "KEGG") |>
    dplyr::select("pathway", "NES") |>
    dplyr::mutate(pathway = stringr::str_replace(.data$pathway, " - Homo.*$", ""))

  y <-
    range(x$NES) |>
    abs() |>
    max()
  lim <- ceiling(y * 100) / 100

  gt::gt(x) |>
    gt::tab_header(
      title = title
    ) |>
    gt::cols_label(
      pathway = "PATHWAY"
    ) |>
    gt::fmt_scientific(
      columns = c("NES")
    ) |>
    gt::data_color(
      columns = .data$NES,
      fn = scales::col_numeric(
        palette =
          grDevices::colorRamp(
            c(clr[1], "white", clr[2]),
            interpolate = "linear"
          ),
        domain = c(-lim, lim)
      )
    ) |>
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = list(
        gt::cells_title(),
        gt::cells_column_labels()
      )
    ) |>
    # gt::tab_style(
    #   style = gt::cell_borders(sides = "bottom"),
    #   locations = list(
    #     gt::cells_body(rows = x$NES == min(x$NES[x$NES > 0]))
    #   )
    # ) |>
    gt::cols_align("center", c(.data$NES)) |>
    gtExtras::gt_theme_538() |>
    gt::opt_table_font(font = "Calibri")
}


# analysis ----------------------------------------------------------------

format_metab_targeted <- function(path) {
  readxl::read_excel(
    path,
    sheet = 2
  ) |>
    dplyr::select(
      id = `Raw File Name`,
      sample = `Sample ID`,
      metabolite = `Compound Name`,
      mz = `Detected Mass`,
      rt = RT,
      area = `Peak Area`
    ) |>
    dplyr::filter(!all(area == 0), .by = metabolite) |>
    dplyr::left_join(wmo::hmdb_mappings, by = "metabolite") |>
    dplyr::mutate(
      metabolite = dplyr::case_when(
        metabolite == "glyceraldehyde 3-phosphate" ~ "GAP",
        metabolite == "glucose 1-phosphate" ~ "G1P",
        metabolite == "glucose 6-phosphate" ~ "G6P",
        metabolite == "sedoheptulose 7-phosphate" ~ "S7P",
        metabolite == "glycerol 3-phosphate" ~ "G3P",
        TRUE ~ metabolite
      )
    ) |>
    dplyr::filter(!is.na(hmdb)) |>
    dplyr::filter(!stringr::str_detect(hmdb, "ISTD")) |>
    dplyr::mutate(
      type = dplyr::case_when(
        sample %in% c("water", "blank") ~ "blank",
        stringr::str_detect(sample, "qc|QC") ~ "qc",
        TRUE ~ "sample"
      ),
      mz = replace(mz, mz == "N/F", NA_real_),
      mz = as.numeric(mz),
      dplyr::across(c(rt, area), \(x) replace(x, x == 0, NA_real_))
    )
}

make_mcti_samples <- function() {
  sample <- as.character(1:60)
  condition <- rep(c("tgfb", "control"), each = 5)
  treatment <- rep(c("DMSO", "VB124", "AZD3965", "AZD/VB", "none"), 2)
  replicate <- rep(1:6, each = 10)

  tibble::tibble(
    sample = sample,
    condition = rep(condition, 6),
    treatment = rep(treatment, 6),
    replicate = as.factor(replicate)
  ) |>
    dplyr::mutate(
      group = stringr::str_c(.data$condition, .data$treatment, sep = "\n")
    ) |>
    refactor() |>
    dplyr::mutate(
      dplyr::across(c("condition", "treatment", "group"), forcats::fct_drop)
    )
}

format_metab_mcti_tar <- function(df, samples) {
  dplyr::left_join(df, samples, by = "sample") |>
    dplyr::select(-"sample") |>
    new_tbl_se(
      a_data = "area",
      f_names = "hmdb",
      f_data = c("metabolite", "mz", "rt"),
      s_names = "id",
      s_data = c("type", "condition", "treatment", "replicate", "group")
    ) |>
    tbl_to_se()
}

make_bleo_samples <- function(animals) {
  sample <- as.character(18:46)
  start_date <- lubridate::as_date(rep("2022-04-01", 29))
  mouse_no <- 1:29

  tibble::tibble(
    sample = sample,
    start_date = start_date,
    mouse_no = mouse_no
  ) |>
    dplyr::left_join(animals, by = c("start_date", "mouse_no")) |>
    dplyr::mutate(
      group = factor(
        group,
        levels = c("Veh", "Bleo", "AZD", "VB"),
        labels = c("Ctl", "Bleo", "AZD", "VB")
      )
    )
}

format_metab_bleo_tar <- function(df, samples) {
  dplyr::left_join(df, samples, by = "sample") |>
    new_tbl_se(
      a_data = "area",
      f_names = "hmdb",
      f_data = c("metabolite", "mz", "rt"),
      s_names = "id",
      s_data = c("type", "idx", "mouse_no", "age", "sex", "group")
    ) |>
    tbl_to_se()
}

plot_bleo_lactate <- function(x) {
  df <-
    purrr::map(x, se_to_tibble) |>
    dplyr::bind_rows(.id = "source") |>
    dplyr::filter(metabolite == "lactate") |>
    dplyr::select(source, group, idx, area) |>
    tidyr::pivot_wider(
      id_cols = c("idx", "group"),
      names_from = "source",
      values_from = "area"
    ) |>
    dplyr::mutate(ratio = lung / plasma) |>
    tidyr::pivot_longer(
      cols = c("plasma", "lung", "ratio"),
      names_to = "measurement",
      values_to = "value"
    ) |>
    dplyr::group_by(measurement) |>
    dplyr::mutate(
      value = value / mean(value[group == "Ctl"]),
      measurement = factor(measurement, levels = c("plasma", "lung", "ratio"))
    )

  stats <-
    df |>
    dplyr::group_by(measurement) |>
    tidyr::nest() |>
    dplyr::mutate(stats = purrr::map(data, \(x) {
      stats::lm(log(value) ~ group, data = x) |>
        emmeans::emmeans(~ group) |>
        emmeans::contrast(
          method = list(
            Ctl = c(-1, 1, 0, 0),
            AZD = c(0, -1, 1, 0),
            VB =  c(0, -1, 0, 1)
          ),
          adjust = "dunnettx"
        ) |>
        broom::tidy() |>
        dplyr::rename(group = contrast) |>
        dplyr::mutate(
          group = factor(group, levels = c("Ctl", "Bleo", "AZD", "VB")),
          label = annot_p(adj.p.value),
          y = Inf,
          vjust = 1.2
        )})) |>
    dplyr::select(measurement, stats) |>
    tidyr::unnest(c(stats)) |>
    dplyr::mutate(measurement = factor(measurement, levels = c("plasma", "lung", "ratio")))

  ggplot2::ggplot(df) +
    ggplot2::facet_wrap(
      ggplot2::vars(measurement),
      # scales = "free_y",
      labeller = ggplot2::labeller(.default = toupper)
    ) +
    ggplot2::aes(
      x = group,
      y = value,
      color = group,
      fill = group
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(color = .data$group),
      stat = "summary",
      fun.data = ggplot2::mean_se,
      width = 0.1,
      linewidth = 0.25,
      position = ggplot2::position_nudge(x = 0.25),
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = .data$group),
      stat = "summary",
      fun = "mean",
      pch = 21,
      size = 1.5,
      color = "white",
      position = ggplot2::position_nudge(x = 0.25),
      show.legend = FALSE
    ) +
    ggbeeswarm::geom_quasirandom(
      pch = 1,
      color = "black",
      width = 0.05,
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = stats,
      ggplot2::aes(
        y = y,
        vjust = vjust,
        label = label
      ),
      color = "black"
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Lactate"
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 6, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      limits = nice_limits,
      label = scales::label_number(scale_cut = scales::cut_short_scale())
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}
