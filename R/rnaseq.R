# rnaseq.R

count_rnaseq <- function() {
  dds <-
    DESeq2::DESeqDataSet(
      rnaseq.lf.tgfb.mcti::se,
      design = ~ replicate + condition * treatment
    ) |>
    DESeq2::DESeq()
  rownames(dds) <- stringr::str_replace_all(rownames(dds), "\\.\\d+", "")
  keep <- rowSums(DESeq2::counts(dds, normalized = TRUE) >= 0.5) >= 4
  dds[keep, ]
}

vst_rnaseq <- function(dds) {
  vsd <- DESeq2::vst(dds, blind = FALSE)
  SummarizedExperiment::assay(vsd) <-
    limma::removeBatchEffect(SummarizedExperiment::assay(vsd), vsd$replicate)

  DESeq2::plotPCA(
    vsd,
    intgroup = c("condition", "treatment"),
    returnData = TRUE
  ) |>
    dplyr::mutate(
      label = stringr::str_replace(group, ":", "\n")
    )
}

plot_rnaseq_pca <- function(pca_data) {
  percent_variance <- round(100 * attr(pca_data, "percentVar"))

  ggplot2::ggplot(pca_data) +
    ggplot2::aes(
      x = PC1,
      y = PC2,
      color = label
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
    ggplot2::geom_point(
      ggplot2::aes(
        fill = label,
      ),
      pch = 21,
      color = "white",
      size = 1.5,
      show.legend = TRUE
    ) +
    ggplot2::labs(
      x = paste0("PC1: ", percent_variance[1], "% variance"),
      y = paste0("PC2: ", percent_variance[2], "% variance")
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      breaks = c("TGFβ\nVeh", "TGFβ\nAZD", "TGFβ\nVB", "TGFβ\nAZD/VB", "Veh", "AZD", "VB", "AZD/VB"),
      labels = c("Veh", "AZD", "VB", "AZD/VB", "Veh", "AZD", "VB", "AZD/VB"),
      values = clrs,
      aesthetics = c("color", "fill"),
      limits = force
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::breaks_extended(n = 5, only.loose = TRUE),
      expand = ggplot2::expansion(c(0.1, 0.1)),
      limits = nice_limits
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 4, only.loose = TRUE),
      expand = ggplot2::expansion(c(0.1, 0.1)),
      limits = nice_limits
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::coord_fixed(clip = "off") +
    ggplot2::theme(
      legend.position = "bottom",
      legend.key.size = ggplot2::unit(1, units = "lines"),
      legend.margin = ggplot2::margin(t = -60),
      legend.text.align = 0
    )
}

plot_counts <- function(dds, symbols) {
  rows <- SummarizedExperiment::rowData(dds)$symbol %in% symbols
  dds[rows, ] |>
    se_to_tibble() |>
    dplyr::select("feature", "symbol", "condition", "treatment", "area") |>
    plot_two_factor(
      x = "treatment",
      y = "area",
      ytitle = "Normalized counts",
      wrap = "symbol"
    )
}

rnaseq_results <- function(dds, con) {
  alpha <- 0.05
  fc <- 1

  annots <-
    tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "row") |>
    dplyr::select(row, entrezid, symbol)

  res <-
    DESeq2::results(
      dds,
      contrast = con,
      alpha = alpha,
      lfcThreshold = log(fc, base = 2),
      parallel = FALSE
    )

  DESeq2::lfcShrink(
    dds,
    res = res,
    type = "ashr"
  ) |>
    tibble::as_tibble(rownames = "row") |>
    dplyr::mutate(stat = res$stat) |>
    dplyr::left_join(annots, by = "row") |>
    dplyr::relocate(symbol, entrezid, .after = "row") |>
    dplyr::arrange(padj)
}

get_msigdb_pathways <- function(species = 'Homo sapiens', category = NULL, subcategory = NULL){
  x <-
    msigdbr::msigdbr(species = species, category = category, subcategory = subcategory) |>
    dplyr::select(gs_name, ensembl_gene) |>
    dplyr::group_by(gs_name) |>
    dplyr::mutate(
      gs_name = ifelse(
        stringr::str_detect(gs_name, 'TARGET_GENES'),
        stringr::str_c(
          'GTRD_',
          stringr::str_replace(gs_name, '_TARGET_GENES', '')
        ),
        gs_name
      )
    ) |>
    tidyr::nest()

  purrr::map(x$data, ~unlist(unname(as.list(.x)), recursive = FALSE)) |>
    rlang::set_names(x$gs_name)
}

plot_vol_rnaseq <- function(
    results,
    gois = NULL,
    title = NULL,
    bins = c(100, 100)
) {
  x <-
    dplyr::filter(results, !is.na(.data$padj)) |>
    dplyr::filter(abs(log2FoldChange) < 10)

  left <-
    x |>
    dplyr::filter(log2FoldChange < 0 & padj < 0.05) |>
    dplyr::filter(!stringr::str_detect(symbol, "LOC\\d+")) |>
    dplyr::slice_min(stat, n = 15)

  right <-
    x |>
    dplyr::filter(log2FoldChange > 0 & padj < 0.05) |>
    dplyr::filter(!stringr::str_detect(symbol, "LOC\\d+")) |>
    dplyr::slice_max(stat, n = 15)

  nudge <- max(0.1 * abs(c(min(x$log2FoldChange), max(x$log2FoldChange))))
  left_nudge <- floor(min(x$log2FoldChange) - nudge)
  right_nudge <- ceiling(max(x$log2FoldChange) + nudge)

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      xlim = c(-Inf, Inf),
      nudge_x = left_nudge - left$log2FoldChange,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = symbol,
        color = symbol %in% gois
      ),
      size = 5/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      xlim = c(-Inf, Inf),
      nudge_x = right_nudge - right$log2FoldChange,
      hjust = 0,
      direction = "y",
      family = "Calibri",
      segment.color = "black",
      show.legend = FALSE
    ) +
    ggplot2::geom_hex(
      bins = bins,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_viridis_c(trans = "log10") +
    ggplot2::scale_color_manual(values = c("black", "darkred")) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log(),
      breaks = scales::breaks_log(4)
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 6),
      expand = ggplot2::expansion(mult = 0.2),
      labels = scales::label_math(2^.x)
    ) +
    ggplot2::labs(
      x = "Fold change",
      y = "Adjusted P value",
      title = title
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}

run_gsea <- function(results, pathways) {
  rnks <-
    results |>
    dplyr::select(row, stat) |>
    dplyr::arrange(stat) |>
    tibble::deframe()

  fgsea::fgsea(
    pathways = pathways,
    stats = rnks,
    nPermSimple = 10000,
    eps = 0
  ) |>
    tibble::as_tibble() |>
    dplyr::filter(.data$padj < 0.05) |>
    dplyr::arrange(desc(NES)) |>
    tidyr::separate(pathway, c("source", "pathway"), "_", extra = "merge")
}

plot_gsea_table <- function(df, title, clr, filename) {
  clr <- clrs[clr]
  x <-
    df |>
    dplyr::select("pathway", "NES") |>
    dplyr::mutate(pathway = stringr::str_replace(.data$pathway, " - Homo.*$", "")) |>
    dplyr::mutate(
      pathway = stringr::str_replace_all(pathway, "_", " "),
      pathway = toupper(pathway)
    )

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
    gt::cols_align("center", c(.data$NES)) |>
    gtExtras::gt_theme_538() |>
    gt::opt_table_font(font = "Calibri")
}

plot_gsea_dot <- function(x) {
  df <-
    dplyr::bind_rows(x, .id = "comparison") |>
    dplyr::filter(padj < 0.05) |>
    dplyr::mutate(
      pathway = stringr::str_replace_all(pathway, "_", " "),
      comparison = factor(
        comparison,
        levels = c("TGFβ\n", "TGFβ +\nAZD", "TGFβ +\nVB", "TGFβ +\nAZD/VB")
      )
    )

  ord <-
    df |>
    dplyr::select(pathway, comparison, NES) |>
    tidyr::pivot_wider(names_from = "comparison", values_from = "NES") |>
    dplyr::arrange(
      dplyr::desc("TGFβ\n"),
      dplyr::desc("TGFβ +\nAZD"),
      dplyr::desc("TGFβ +\nVB"),
      dplyr::desc("TGFβ +\nAZD/VB")
    ) |>
    dplyr::pull(pathway)

  df |>
    dplyr::mutate(pathway = factor(pathway, levels = ord)) |>
    ggplot2::ggplot() +
    ggplot2::aes(
      x = comparison,
      y = pathway
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = NES
      ),
      color = "white",
      size = 2,
      pch = 21
    ) +
    ggplot2::scale_fill_distiller(
      type = "div",
      palette = "RdBu",
      limits = c(-2, 2),
      oob = scales::squish
    ) +
    ggplot2::scale_y_discrete(
      limits = rev
    ) +
    ggplot2::scale_x_discrete(
      position = "top"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL
    ) +
    ggplot2::guides(
      size = "none",
      fill = ggplot2::guide_colorbar(barheight = 5)
    ) +
    theme_plot() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(hjust = 0),
      axis.text.x = ggplot2::element_text(vjust = 0),
      panel.grid.major.y = ggplot2::element_line(
        linewidth = ggplot2::rel(0.25),
        color = "grey80"
      ),
      legend.position = "right",
      legend.text.align = 1
    ) +
    NULL
}

plot_edge <- function(x, dds) {
  targets <-
    dplyr::bind_rows(x, .id = "treatment") |>
    dplyr::filter(pathway == "EPITHELIAL_MESENCHYMAL_TRANSITION") |>
    dplyr::select(treatment, leadingEdge) |>
    tidyr::unnest(c(leadingEdge)) |>
    dplyr::group_by(leadingEdge) |>
    dplyr::count() |>
    dplyr::filter(n == 3) |>
    dplyr::pull(leadingEdge)

  x <-
    dds[targets, ] |>
    DESeq2::counts(normalize = TRUE) |>
    tibble::as_tibble(rownames = "row") |>
    tidyr::pivot_longer(
      cols = !row,
      names_to = "sample",
      values_to = "count"
    ) |>
    dplyr::left_join(
      tibble::as_tibble(SummarizedExperiment::rowData(dds), rownames = "row"),
      by = "row",
      copy = TRUE
    ) |>
    dplyr::left_join(
      SummarizedExperiment::colData(dds),
      by = c("sample" = "names"),
      copy = TRUE
    ) |>
    dplyr::select(
      symbol, sample, replicate, condition, treatment, count
    ) |>
    dplyr::mutate(group = stringr::str_c(condition, treatment, sep = "\n")) |>
    refactor() |>
    dplyr::filter(group %in% c("Ctl\nVeh", "TGFβ\nVeh", "TGFβ\nAZD", "TGFβ\nVB"))

  ggplot2::ggplot(x) +
    ggplot2::facet_wrap(
      ggplot2::vars(symbol),
      scales = "free_y",
      nrow = 2
    ) +
    ggplot2::aes(
      x = group,
      y = count
    ) +
    ggplot2::geom_point(
      pch = 1,
      color = "black",
      size = 0.5,
      stroke = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        color = group
      ),
      stat = "summary",
      fun.data = "mean_se",
      width = 0.25,
      linewidth = 0.25,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(
        fill = group
      ),
      size = 1.5,
      stat = "summary",
      fun = "mean",
      pch = 21,
      color = "white"
    ) +
    ggplot2::scale_color_manual(
      values = clrs,
      aesthetics = c("color", "fill")
    ) +
    ggplot2::scale_y_continuous(
      breaks = scales::breaks_extended(n = 5, only.loose = TRUE),
      expand = ggplot2::expansion(c(0, 0)),
      labels = scales::label_number(scale_cut = c(0, K = 1000)),
      limits = nice_limits
    ) +
    ggplot2::labs(
      title = "EMT Leading Edge Genes",
      x = NULL,
      y = "Count",
      color = NULL,
      fill = NULL
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(nrow = 2, byrow = TRUE)) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      legend.position = c(0.87, 0.22),
      legend.key.size = ggplot2::unit(1, "lines")
    )
}

plot_gsea_vol <- function(deg, pathways, pathway) {
  title <-
    stringr::str_replace_all(pathway, "_", " ") |>
    stringr::str_replace("HALLMARK ", "")

  x <-
    deg |>
    dplyr::filter(row %in% pathways[[pathway]])

  left <-
    x |>
    dplyr::filter(log2FoldChange < 0 & padj < 0.05) |>
    dplyr::slice_min(stat, n = 10)

  right <-
    x |>
    dplyr::filter(log2FoldChange > 0 & padj < 0.05) |>
    dplyr::slice_max(stat, n = 10)

  nudge <- max(0.1 * abs(c(min(x$log2FoldChange), max(x$log2FoldChange))))
  left_nudge <- min(x$log2FoldChange) - nudge
  right_nudge <- max(x$log2FoldChange) + nudge

  ggplot2::ggplot(df) +
    ggplot2::aes(
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = symbol
      ),
      color = "black",
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = left_nudge - left$log2FoldChange,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = symbol
      ),
      color = "black",
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      segment.color = "black",
      nudge_x = right_nudge - right$log2FoldChange,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      # ggplot2::aes(fill = factor(sign)),
      pch = 21,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log(),
      breaks = scales::breaks_log(6)
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 7)
    ) +
    ggplot2::labs(
      x = "Fold Change",
      y = "Adjusted P value",
      title = title
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.y = ggplot2::element_text(hjust = 0)
    ) +
    NULL
}

prep_counts <- function(dds) {
  # remove duplicates
  ids <-
    SummarizedExperiment::rowData(dds) |>
    tibble::as_tibble(rownames = "row") |>
    dplyr::filter(!(is.na(symbol) | symbol == "")) |>
    dplyr::group_by(symbol) |>
    dplyr::slice_max(baseMean) |>
    dplyr::select(row, symbol) |>
    tibble::deframe()

  # make matrix
  dds[rownames(dds) %in% names(ids), ] |>
    DESeq2::counts(normalized = TRUE) |>
    {\(x) magrittr::set_rownames(x, ids[rownames(x)])}()
}

plot_progeny <- function(df) {
  palette_length <- 100

  x <-
    df |>
    dplyr::mutate(
      name = .data$condition,
      condition = stringr::str_extract(.data$name, "(Ctl|TGFb)"),
      condition = factor(condition, levels = c("Ctl", "TGFb")),
      treatment = stringr::str_extract(.data$name, "(DMSO|AZD|VB|Dual)"),
      treatment = factor(treatment, levels = c("DMSO", "AZD", "VB", "Dual"))
    ) |>
    dplyr::arrange(source, condition, treatment, name)

  mat <-
    x |>
    tidyr::pivot_wider(
      id_cols = "source",
      names_from = "name",
      values_from = "score"
    ) |>
    tibble::column_to_rownames("source") |>
    as.matrix() |>
    t() |>
    scale() |>
    t()

  annotation_col <-
    x |>
    dplyr::select(name, condition, treatment) |>
    dplyr::distinct() |>
    refactor() |>
    tibble::column_to_rownames("name")

  annotation_colors <-
    list(
      condition = c(
        Ctl = clrs[["Ctl"]],
        TGFβ = clrs[["TGFβ"]]
      ),
      treatment = c(
        Veh = clrs[["TGFβ\nVeh"]],
        AZD = clrs[["TGFβ\nAZD"]],
        VB = clrs[["TGFβ\nVB"]],
        `AZD/VB` = clrs[["TGFβ\nAZD/VB"]]
      )
    )

  pheatmap::pheatmap(
    mat,
    color = colorspace::diverging_hcl(100, "Blue-Red 3"),
    breaks = seq(-3, 3, length.out = palette_length),
    cluster_cols = FALSE,
    border_color = NA,
    cellwidth = 10,
    cellheight = 10,
    legend = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = FALSE
  )
}

# run_tfea <- function(dds){
#   # vst if using scale as method
#   # dds <- DESeq2::vst(dds, blind = FALSE)
#
#   # remove duplicates
#   ids <-
#     SummarizedExperiment::rowData(dds) |>
#     tibble::as_tibble() |>
#     dplyr::filter(!(is.na(symbol) | symbol == "")) |>
#     dplyr::group_by(symbol) |>
#     dplyr::slice_max(baseMean) |>
#     dplyr::select(gene_id, symbol) |>
#     tibble::deframe()
#
#   # make matrix
#   mat <-
#     dds[rownames(dds) %in% names(ids), ] |>
#     SummarizedExperiment::assay() |>
#     {\(x) magrittr::set_rownames(x, ids[rownames(x)])}()
#
#   # run viper
#   regulons <-
#     dorothea::dorothea_hs |>
#     dplyr::filter(confidence %in% c("A", "B"))
#
#   dorothea::run_viper(
#     input = mat,
#     regulons = regulons,
#     options = list(
#       method = "rank",
#       minsize = 1,
#       nes = TRUE,
#       cores = 4,
#       verbose = FALSE
#     ),
#     tidy = FALSE
#   )
# }

run_tfea <- function(mat, net) {
  decoupleR::decouple(mat, net) |>
    dplyr::filter(statistic == "consensus")
}

plot_tf <- function(df) {
  # prevent Rplots.pdf generation
  if(!interactive()) pdf(NULL)

  palette_length <- 100

  x <-
    df |>
    dplyr::mutate(
      name = .data$condition,
      condition = stringr::str_extract(.data$name, "(Ctl|TGFb)"),
      condition = factor(condition, levels = c("Ctl", "TGFb")),
      treatment = stringr::str_extract(.data$name, "(DMSO|AZD|VB|Dual)"),
      treatment = factor(treatment, levels = c("DMSO", "AZD", "VB", "Dual"))
    ) |>
    # dplyr::filter(condition != "Ctl") |>
    dplyr::arrange(source, condition, treatment, name)

  mat <-
    x |>
    tidyr::pivot_wider(
      id_cols = "source",
      names_from = "name",
      values_from = "score"
    ) |>
    tibble::column_to_rownames("source") |>
    as.matrix() |>
    t() |>
    scale() |>
    t()

  tfs <-
    df |>
    dplyr::summarise(std = sd(score), .by = "source") |>
    dplyr::arrange(-abs(std)) |>
    dplyr::slice_head(n = 25) |>
    dplyr::pull(source)

  annotation_col <-
    x |>
    dplyr::select(name, condition, treatment) |>
    dplyr::distinct() |>
    refactor() |>
    tibble::column_to_rownames("name")

  annotation_colors <-
    list(
      condition = c(
        Ctl = clrs[["Ctl"]],
        TGFβ = clrs[["TGFβ"]]
      ),
      treatment = c(
        Veh = clrs[["TGFβ\nVeh"]],
        AZD = clrs[["TGFβ\nAZD"]],
        VB = clrs[["TGFβ\nVB"]],
        `AZD/VB` = clrs[["TGFβ\nAZD/VB"]]
      )
    )

  pheatmap::pheatmap(
    mat[tfs, ],
    color = colorspace::diverging_hcl(100, "Blue-Red 3"),
    breaks = seq(-3.5, 3.5, length.out = palette_length),
    border_color = NA,
    cluster_cols = FALSE,
    legend = TRUE,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    show_colnames = FALSE
  )
}

score_tf <- function(deg, net) {
  mat <-
    deg |>
    dplyr::filter(!is.na(symbol)) |>
    dplyr::select(symbol, stat) |>
    dplyr::group_by(symbol) |>
    dplyr::slice_max(order_by = abs(stat), n = 1) |>
    tibble::column_to_rownames("symbol")

  decoupleR::decouple(
    mat = mat,
    network = net
  ) |>
    dplyr::filter(statistic == "consensus") |>
    dplyr::select(source, score, pval = p_value)
}

plot_vol_tfea <- function(
    tfea,
    fills = c("Ctl", "IPF"),
    title = NULL,
    mois = NULL
) {
  left <-
    tfea |>
    dplyr::filter(score < 0 & pval < 0.1) |>
    dplyr::slice_min(score * -log10(pval), n = 10)

  right <-
    tfea |>
    dplyr::filter(score > 0 & pval < 0.1) |>
    dplyr::slice_max(score * -log10(pval), n = 10)

  nudge <- max(0.2 * abs(c(min(tfea$score), max(tfea$score))))
  left_nudge <- min(tfea$score) - nudge
  right_nudge <- max(tfea$score) + nudge

  ggplot2::ggplot(tfea) +
    ggplot2::aes(
      x = .data$score,
      y = .data$pval
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = source,
        color = source %in% mois
      ),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = left_nudge - left$score,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = source,
        color = source %in% mois
      ),
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      segment.color = "black",
      nudge_x = right_nudge - right$score,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = subset(tfea, pval > 0.1),
      pch = 21,
      color = "white",
      fill = "grey80"
    ) +
    ggplot2::geom_point(
      data = subset(tfea, score > 0 & pval < 0.1),
      pch = 21,
      color = "white",
      fill = clrs[[fills[[2]]]]
    ) +
    ggplot2::geom_point(
      data = subset(tfea, score < 0 & pval < 0.1),
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
      breaks = scales::pretty_breaks(n = 7)
    ) +
    ggplot2::labs(
      x = "Score",
      y = "P value",
      title = title
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.text.y = ggplot2::element_text(hjust = 0)
    ) +
    NULL
}

plot_tf_targets <- function(tf, targets, deg) {
  x <-
    targets |>
    dplyr::filter(source == tf) |>
    dplyr::inner_join(deg, by = c("target" = "symbol")) |>
    dplyr::mutate(sign = sign(mor) * sign(stat))

  left <-
    x |>
    dplyr::filter(sign == 1 & log2FoldChange < 0 & padj < 0.05) |>
    dplyr::slice_min(stat, n = 10)

  right <-
    x |>
    dplyr::filter(sign == 1 & log2FoldChange > 0 & padj < 0.05) |>
    dplyr::slice_max(stat, n = 10)

  nudge <- max(0.1 * abs(c(min(x$log2FoldChange), max(x$log2FoldChange))))
  left_nudge <- min(x$log2FoldChange) - nudge
  right_nudge <- max(x$log2FoldChange) + nudge

  ggplot2::ggplot(x) +
    ggplot2::aes(
      x = log2FoldChange,
      y = padj
    ) +
    ggrepel::geom_text_repel(
      data = left,
      ggplot2::aes(
        label = target
      ),
      color = "black",
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      nudge_x = left_nudge - left$log2FoldChange,
      hjust = 0,
      segment.color = "black",
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggrepel::geom_text_repel(
      data = right,
      ggplot2::aes(
        label = target
      ),
      color = "black",
      size = 6/ggplot2::.pt,
      max.overlaps = 20,
      segment.size = 0.1,
      segment.color = "black",
      nudge_x = right_nudge - right$log2FoldChange,
      hjust = 1,
      direction = "y",
      family = "Calibri",
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = factor(sign)),
      pch = 21,
      color = "white",
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_manual(values = c("grey50", "darkred")) +
    ggplot2::scale_y_continuous(
      trans = c("log10", "reverse"),
      labels = scales::label_log(),
      breaks = scales::breaks_log(4)
    ) +
    ggplot2::scale_x_continuous(
      breaks = scales::pretty_breaks(n = 4),
      # limits = c(left_nudge - 0.1, right_nudge + 0.1),
      labels = scales::label_math(2^.x)
    ) +
    ggplot2::labs(
      x = "Fold change",
      y = "Adjusted P value",
      title = tf
    ) +
    ggplot2::coord_cartesian(clip = "off") +
    theme_plot() +
    NULL
}
