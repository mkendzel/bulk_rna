# Overview, DE and expression figures for all 25 contrasts.
# GSEA figures are in 6_GSEA_figures.R; pathway redundancy QC in 7_pathway_QC.R.
# Run from the repo root.
#
# Figures are driven by contrast_registry. A contrast with nothing significant
# gets an annotated empty panel rather than no file.

# ---- Library ----
library(ggplot2)
library(ggrepel)
library(ggvenn)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(grid)
library(pheatmap)
library(patchwork)
library(openxlsx)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Data load ----
# logcpm is voom's E matrix (log2 counts per million)
tt          <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)), EXPERIMENTS)
logcpm      <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("logcpm_", e), dir = dir_rds)), EXPERIMENTS)
sample_meta <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("sample_meta_", e), dir = dir_rds)), EXPERIMENTS)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)

# ---- Coverage check ----
# Warn when the registry lists a contrast that is absent from tt. 3_DE_limma.R
# drops contrasts whose groups lost all samples but checkpoints the full registry.
coverage <- purrr::map_dfr(EXPERIMENTS, function(e) {
  expect <- contrasts_for(e)
  tibble::tibble(
    experiment = e,
    n_expected = length(expect),
    n_toptable = sum(expect %in% names(tt[[e]])),
    missing    = paste(setdiff(expect, names(tt[[e]])), collapse = ", ")
  )
})
print(coverage)

if (any(coverage$n_toptable != coverage$n_expected)) {
  warning("Contrasts missing from tt - re-run 3_DE_limma.R:\n",
          paste(coverage$missing[nzchar(coverage$missing)], collapse = "\n"),
          call. = FALSE)
}

# Contrasts present in tt for one experiment, in registry order.
contrasts_present <- function(e) intersect(contrasts_for(e), names(tt[[e]]))

# All toptables of one experiment stacked, with registry metadata joined and the
# Up/Down/NS call applied. Feeds the volcano grid, MA plots, p-value histograms
# and the DE count barplot.
long_tt <- function(e) {
  purrr::map_dfr(contrasts_present(e), function(cn) {
    dplyr::mutate(tt[[e]][[cn]], contrast = cn)
  }) |>
    dplyr::left_join(
      dplyr::select(contrast_registry, name, label, type, min_n),
      by = c("contrast" = "name")
    ) |>
    dplyr::mutate(
      label     = factor(label, levels = contrast_labels(e)),
      direction = de_direction(logFC, adj.P.Val)
    )
}

# Facet grid geometry scales with contrast count: 23 spinal, 2 cort.
facet_dims <- function(n, ncol = 4, w_per = 4, h_per = 3.2) {
  ncol <- min(ncol, n)
  list(ncol = ncol, width = max(7, ncol * w_per),
       height = max(5, ceiling(n / ncol) * h_per))
}

# =============================================================================
# Overview
# =============================================================================
dir_overview <- file.path("overview")

# ---- PCA ----
# prcomp once per experiment, shared by the PC1/2, PC3/4 and scree plots.
# Non-finite and zero-variance genes are dropped before ranking, so the top-n
# selection cannot pick up an NA.
pca_fit <- function(lc, n_top = 500) {
  v    <- matrixStats::rowVars(lc)
  names(v) <- rownames(lc)
  v    <- v[is.finite(v) & v > 0]
  keep <- utils::head(names(sort(v, decreasing = TRUE)), min(n_top, length(v)))
  p    <- stats::prcomp(t(lc[keep, , drop = FALSE]), scale. = FALSE)
  list(x = p$x, pct = round(100 * p$sdev^2 / sum(p$sdev^2), 1), n_used = length(keep))
}

pca_plot <- function(fit, sm, pcs = c(1, 2), title = NULL) {

  if (ncol(fit$x) < max(pcs)) {
    return(empty_plot(title, note = paste0("Fewer than PC", max(pcs), " available")))
  }

  d <- tibble::tibble(
    sample = rownames(fit$x),
    xx     = fit$x[, pcs[1]],
    yy     = fit$x[, pcs[2]]
  ) |>
    dplyr::left_join(dplyr::select(sm, sample, line, stim, condition), by = "sample")

  centroids <- d |>
    dplyr::group_by(condition, line, stim) |>
    dplyr::summarise(cx = mean(xx), cy = mean(yy), .groups = "drop")

  d <- dplyr::left_join(d, centroids, by = c("condition", "line", "stim"))

  ggplot(d, aes(xx, yy)) +
    # Legs from each sample to its group centroid
    geom_segment(aes(xend = cx, yend = cy, colour = condition),
                 alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
    geom_point(aes(fill = condition, shape = stim),
               size = 3.2, colour = "black", stroke = 0.5) +
    geom_point(data = centroids, aes(cx, cy, fill = condition),
               shape = 21, size = 5, colour = "black", stroke = 1.1,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = cond_cols, drop = TRUE) +
    scale_colour_manual(values = cond_cols, drop = TRUE) +
    scale_shape_manual(values = stim_shapes, drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = paste0("PC", pcs[1], ": ", fit$pct[pcs[1]], "% variance"),
         y = paste0("PC", pcs[2], ": ", fit$pct[pcs[2]], "% variance"),
         title = title) +
    theme_classic()
}

scree_plot <- function(fit, n_pc = 10, title = NULL) {
  n <- min(n_pc, length(fit$pct))
  d <- tibble::tibble(
    pc  = factor(paste0("PC", seq_len(n)), levels = paste0("PC", seq_len(n))),
    pct = fit$pct[seq_len(n)],
    cum = cumsum(fit$pct)[seq_len(n)]
  )

  ggplot(d, aes(pc, pct)) +
    geom_col(fill = "grey70", colour = "black", width = 0.7) +
    geom_line(aes(y = cum, group = 1), colour = "#D62728", linewidth = 0.7) +
    geom_point(aes(y = cum), colour = "#D62728", size = 2) +
    geom_text(aes(label = sprintf("%.0f", pct)), vjust = -0.4, size = 3) +
    labs(x = NULL, y = "% variance (bars)   |   cumulative (line)", title = title) +
    theme_classic()
}

for (e in EXPERIMENTS) {

  fit <- pca_fit(logcpm[[e]])

  save_fig(pca_plot(fit, sample_meta[[e]], c(1, 2),
                    sprintf("PCA - %s (top %d variable genes)", e, fit$n_used)),
           paste0("pca_", e), width = 8, height = 6, subdir = dir_overview)

  save_fig(pca_plot(fit, sample_meta[[e]], c(3, 4),
                    sprintf("PCA PC3/PC4 - %s (top %d variable genes)", e, fit$n_used)),
           paste0("pca_pc34_", e), width = 8, height = 6, subdir = dir_overview)

  save_fig(scree_plot(fit, 10, paste0("Scree - ", e)),
           paste0("pca_scree_", e), width = 7, height = 5, subdir = dir_overview)
}

# ---- Sample-sample correlation ----
# Pearson over every filtered gene, per experiment, on the voom matrix the DE
# used. Clustered, unlike figures/qc/qc_sample_clustering.png which is RNAseqQC
# on the vst of all 40 samples pooled.
for (e in EXPERIMENTS) {

  lc <- logcpm[[e]]
  sm <- sample_meta[[e]]
  lc <- lc[matrixStats::rowSds(lc) > 0, , drop = FALSE]

  cm <- stats::cor(lc, method = "pearson")

  if (!all(is.finite(cm))) {
    message("Non-finite correlations for ", e, " - skipping sample correlation heatmap")
    next
  }

  ann <- data.frame(line = sm$line, stim = sm$stim, row.names = sm$sample)
  ann_colors <- list(
    line = line_cols[levels(droplevels(sm$line))],
    stim = setNames(grDevices::grey.colors(nlevels(droplevels(sm$stim)),
                                           start = 0.2, end = 0.9),
                    levels(droplevels(sm$stim)))
  )

  ph <- pheatmap::pheatmap(
    cm, silent = TRUE,
    clustering_distance_rows = stats::as.dist(1 - cm),
    clustering_distance_cols = stats::as.dist(1 - cm),
    annotation_col    = ann,
    annotation_colors = ann_colors,
    display_numbers   = FALSE,
    main = sprintf("Sample correlation (Pearson, %d genes) - %s", nrow(lc), e)
  )

  save_fig(ph, paste0("sample_cor_", e), width = 9, height = 8, subdir = dir_overview)
}

# ---- Heatmap: most variable genes ----
for (e in EXPERIMENTS) {

  lc  <- logcpm[[e]]
  sm  <- sample_meta[[e]]
  map <- symbol_map(tt[[e]])

  v <- matrixStats::rowVars(lc)
  names(v) <- rownames(lc)
  v <- v[is.finite(v) & v > 0]
  top_ens <- utils::head(names(sort(v, decreasing = TRUE)), min(50, length(v)))

  mat <- lc[top_ens, , drop = FALSE]
  rownames(mat) <- label_genes(top_ens, map)

  mat_scaled <- zscore_rows(mat)

  if (!heatmap_ready(mat_scaled)) {
    message("Too few variable genes for ", e, " - skipping top-variable heatmap")
    next
  }

  desired_order <- sm |> dplyr::arrange(condition, rep) |> dplyr::pull(sample)
  mat_scaled    <- mat_scaled[, desired_order, drop = FALSE]

  ann <- data.frame(line = sm$line, stim = sm$stim, row.names = sm$sample)
  ann <- ann[desired_order, , drop = FALSE]

  ann_colors <- list(
    line = line_cols[levels(droplevels(sm$line))],
    stim = setNames(grDevices::grey.colors(nlevels(droplevels(sm$stim)),
                                           start = 0.2, end = 0.9),
                    levels(droplevels(sm$stim)))
  )

  ph <- pheatmap::pheatmap(
    mat_scaled, silent = TRUE,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = ann,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    main = paste0("Top ", nrow(mat_scaled), " variable genes - ", e)
  )

  save_fig(ph, paste0("heatmap_top50_variable_", e),
           width = 10, height = 9, subdir = dir_overview)
}

# ---- log2FC heatmap across contrasts ----
top_n_genes <- 250

for (e in EXPERIMENTS) {

  cns     <- contrasts_present(e)
  res_use <- tt[[e]][cns]

  sig_tbl <- purrr::map_dfr(res_use, function(df) {
    df[is_de(df$logFC, df$adj.P.Val), c("ensembl_id", "logFC")]
  })

  if (nrow(sig_tbl) == 0) {
    save_fig(empty_plot(paste0("Top DE genes, log2 fold change - ", e),
                        note = sprintf("No gene passes adj.P.Val < %s and |logFC| > %s",
                                       P_CUT, LFC_CUT)),
             paste0("heatmap_log2FC_", e), width = 8, height = 5, subdir = dir_overview)
    next
  }

  top_genes <- sig_tbl |>
    dplyr::group_by(ensembl_id) |>
    dplyr::summarise(max_abs_lfc = max(abs(logFC)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_abs_lfc)) |>
    dplyr::slice_head(n = top_n_genes) |>   # clamps to available rows
    dplyr::pull(ensembl_id)

  # Non-significant cells become NA and are drawn black
  lfc_mat <- vapply(res_use, function(df) {
    idx <- match(top_genes, df$ensembl_id)
    lfc <- df$logFC[idx]
    lfc[!is_de(lfc, df$adj.P.Val[idx])] <- NA_real_
    lfc
  }, numeric(length(top_genes)))

  reg_e <- contrast_registry |>
    dplyr::filter(name %in% cns) |>
    dplyr::arrange(type, line, ref_line, stim)

  lfc_mat <- lfc_mat[, reg_e$name, drop = FALSE]
  colnames(lfc_mat) <- reg_e$label
  rownames(lfc_mat) <- label_genes(top_genes, symbol_map(res_use))

  ann <- data.frame(type = reg_e$type, row.names = reg_e$label)

  # NA -> 0 for the clustering distance only
  lfc_for_cluster <- lfc_mat
  lfc_for_cluster[is.na(lfc_for_cluster)] <- 0

  ph <- pheatmap::pheatmap(
    lfc_mat, silent = TRUE,
    cluster_rows      = if (heatmap_ready(lfc_for_cluster))
                          stats::hclust(stats::dist(lfc_for_cluster)) else FALSE,
    cluster_cols      = FALSE,
    annotation_col    = ann,
    annotation_colors = list(type = type_cols[levels(droplevels(reg_e$type))]),
    na_col            = "black",
    show_colnames     = TRUE,
    show_rownames     = FALSE,
    angle_col         = 45,
    main = paste0("Top ", nrow(lfc_mat), " DE genes, log2 fold change - ", e,
                  "  (black = not significant)")
  )

  save_fig(ph, paste0("heatmap_log2FC_", e),
           width = max(8, 3 + nrow(reg_e) * 0.42), height = 10, subdir = dir_overview)
}

# ---- DE gene counts across contrasts ----
# Recomputed from tt rather than read back from results/de_gene_counts.tsv,
# which would coerce the registry factors to character.
for (e in EXPERIMENTS) {

  # Built over the full registry so a contrast missing from tt still gets a
  # labelled slot. scale_x_discrete(drop = FALSE) does not work here: with
  # facet_grid(scales = "free_x") every facet would list every contrast.
  counts <- purrr::map_dfr(contrasts_for(e), function(cn) {
    if (!cn %in% names(tt[[e]])) {
      return(tibble::tibble(contrast = cn, up = NA_integer_, down = NA_integer_))
    }
    d <- de_direction(tt[[e]][[cn]]$logFC, tt[[e]][[cn]]$adj.P.Val)
    tibble::tibble(contrast = cn, up = sum(d == "Up"), down = sum(d == "Down"))
  }) |>
    dplyr::left_join(dplyr::select(contrast_registry, name, label, type, min_n),
                     by = c("contrast" = "name")) |>
    dplyr::mutate(label = factor(label, levels = contrast_labels(e)))

  plot_df <- counts |>
    tidyr::pivot_longer(c(up, down), names_to = "direction", values_to = "n") |>
    dplyr::mutate(
      signed    = ifelse(direction == "up", n, -n),
      direction = factor(ifelse(direction == "up", "Up", "Down"), levels = DE_LEVELS)
    )

  p <- ggplot(plot_df, aes(label, signed, fill = direction)) +
    geom_col(colour = "black", width = 0.72, linewidth = 0.25) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_text(aes(label = ifelse(is.na(n), "", ifelse(n == 0, "0", format(n, big.mark = ","))),
                  vjust = ifelse(signed >= 0, -0.4, 1.3)),
              size = 2.5, na.rm = TRUE) +
    # Unreplicated contrasts flagged once, above the panel
    geom_text(data = dplyr::distinct(counts, label, type, min_n),
              aes(x = label, y = Inf, label = ifelse(min_n == 1, "n = 1", "")),
              inherit.aes = FALSE, vjust = 1.4, size = 2.6, colour = "#B22222") +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    scale_y_continuous(labels = function(z) format(abs(z), big.mark = ",")) +
    labs(x = NULL, y = "DE genes (up above, down below)",
         title = paste0("Differentially expressed genes - ", e),
         subtitle = sprintf("adj.P.Val < %s and |log2FC| > %s", P_CUT, LFC_CUT)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          strip.background = element_rect(fill = "grey92", colour = NA))

  save_fig(p, paste0("de_counts_", e),
           width = max(7, 2.2 + nrow(counts) * 0.62), height = 6, subdir = dir_overview)
}

# =============================================================================
# Per-contrast DE figures
# =============================================================================

# ---- Volcano ----
volcano_plot <- function(df, title = NULL, subtitle = NULL, caption = NULL,
                         label_n = 0) {

  d <- dplyr::mutate(df,
    neglog10_padj = -log10(adj.P.Val),
    direction     = de_direction(logFC, adj.P.Val)
  )

  n_up   <- sum(d$direction == "Up",   na.rm = TRUE)
  n_down <- sum(d$direction == "Down", na.rm = TRUE)

  p <- ggplot(d, aes(logFC, neglog10_padj)) +
    geom_point(aes(colour = direction), alpha = 0.6, size = 1) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
    geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey40") +
    # Inf placement keeps the counts on-panel when the y axis is near-degenerate
    annotate("text", x = Inf, y = Inf, label = paste0("Up: ", n_up),
             hjust = 1.1, vjust = 1.5, colour = DE_COLS[["Up"]], size = 3.2) +
    annotate("text", x = -Inf, y = Inf, label = paste0("Down: ", n_down),
             hjust = -0.1, vjust = 1.5, colour = DE_COLS[["Down"]], size = 3.2) +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = "log2 fold change", y = "-log10 adjusted p-value") +
    theme_classic() +
    theme(plot.caption = element_text(colour = "#B22222", hjust = 0))

  if (label_n > 0) {
    lab <- d |>
      dplyr::filter(direction != "NS", !is.na(SYMBOL), SYMBOL != "") |>
      dplyr::slice_min(adj.P.Val, n = label_n, with_ties = FALSE)
    if (nrow(lab) > 0) {
      p <- p + ggrepel::geom_text_repel(
        data = lab, aes(label = SYMBOL), size = 2.6, fontface = "italic",
        max.overlaps = 20, min.segment.length = 0, segment.size = 0.3,
        segment.colour = "grey60")
    }
  }
  p
}

for (e in EXPERIMENTS) {

  sub_v <- file.path("de", "volcano", e)

  for (cn in contrasts_present(e)) {
    r <- reg_of(cn)
    save_fig(
      volcano_plot(tt[[e]][[cn]], title = r$label, subtitle = cn,
                   caption = if (isTRUE(r$min_n == 1)) "Unreplicated (n = 1) - exploratory" else NULL,
                   label_n = 15),
      paste0("volcano_", cn), width = 7, height = 6, subdir = sub_v)
  }

  # Comparative grid
  d  <- long_tt(e)
  fd <- facet_dims(nlevels(droplevels(d$label)))

  p <- ggplot(d, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(colour = direction), alpha = 0.5, size = 0.5) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey55") +
    geom_hline(yintercept = -log10(P_CUT), linetype = "dashed", colour = "grey55") +
    facet_wrap(~ label, ncol = fd$ncol, drop = FALSE) +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = paste0("Volcano - ", e), x = "log2 fold change",
         y = "-log10 adjusted p-value") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5))

  save_fig(p, paste0("volcano_grid_", e), width = fd$width, height = fd$height,
           subdir = file.path("de", "volcano"))
}

# ---- MA plots ----
ma_plot <- function(df, title = NULL, subtitle = NULL, caption = NULL) {
  d <- dplyr::mutate(df, direction = de_direction(logFC, adj.P.Val))

  ggplot(d, aes(AveExpr, logFC)) +
    geom_point(aes(colour = direction), alpha = 0.6, size = 1) +
    geom_hline(yintercept = 0, linewidth = 0.4) +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey40") +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = title, subtitle = subtitle, caption = caption,
         x = "Average expression (log2 CPM)", y = "log2 fold change") +
    theme_classic() +
    theme(plot.caption = element_text(colour = "#B22222", hjust = 0))
}

for (e in EXPERIMENTS) {

  sub_m <- file.path("de", "ma", e)

  for (cn in contrasts_present(e)) {
    r <- reg_of(cn)
    save_fig(
      ma_plot(tt[[e]][[cn]], title = r$label, subtitle = cn,
              caption = if (isTRUE(r$min_n == 1)) "Unreplicated (n = 1) - exploratory" else NULL),
      paste0("ma_", cn), width = 7, height = 6, subdir = sub_m)
  }

  d  <- long_tt(e)
  fd <- facet_dims(nlevels(droplevels(d$label)))

  p <- ggplot(d, aes(AveExpr, logFC)) +
    geom_point(aes(colour = direction), alpha = 0.5, size = 0.5) +
    geom_hline(yintercept = 0, linewidth = 0.35) +
    geom_hline(yintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", colour = "grey55") +
    facet_wrap(~ label, ncol = fd$ncol, drop = FALSE) +
    scale_colour_manual(values = DE_COLS, drop = FALSE, name = NULL) +
    labs(title = paste0("MA - ", e), x = "Average expression (log2 CPM)",
         y = "log2 fold change") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5))

  save_fig(p, paste0("ma_grid_", e), width = fd$width, height = fd$height,
           subdir = file.path("de", "ma"))
}

# ---- p-value histograms ----
# One facet per contrast, annotated with the observed/expected gene count at
# P_CUT (ratio) and the estimated null fraction (pi0).
for (e in EXPERIMENTS) {

  d <- long_tt(e)

  stats_df <- d |>
    dplyr::group_by(label) |>
    dplyr::summarise(
      n_genes = sum(!is.na(P.Value)),
      ratio   = sum(P.Value < P_CUT, na.rm = TRUE) / (P_CUT * n_genes),
      pi0     = min(1, 2 * mean(P.Value > 0.5, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(txt = sprintf("ratio %.2f   pi0 %.2f", ratio, pi0))

  fd <- facet_dims(nlevels(droplevels(d$label)))

  p <- ggplot(d, aes(P.Value)) +
    geom_histogram(binwidth = 0.025, boundary = 0,
                   fill = "grey70", colour = "black", linewidth = 0.15) +
    # Null-uniform expectation per bin
    geom_hline(data = stats_df, aes(yintercept = n_genes * 0.025),
               linetype = "dashed", colour = "#D62728", linewidth = 0.4) +
    geom_text(data = stats_df, aes(x = 0.98, y = Inf, label = txt),
              inherit.aes = FALSE, hjust = 1, vjust = 1.6, size = 2.6, colour = "grey20") +
    facet_wrap(~ label, ncol = fd$ncol, scales = "free_y", drop = FALSE) +
    labs(title = paste0("Raw p-value distribution - ", e),
         subtitle = paste0("dashed line = uniform null;  ",
                           "ratio = observed / expected genes at p < ", P_CUT,
                           ";  pi0 = estimated null fraction"),
         x = "Raw p-value", y = "Genes") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          strip.text = element_text(size = 7.5))

  save_fig(p, paste0("pvalue_hist_", e), width = fd$width, height = fd$height,
           subdir = file.path("de", "diagnostics"))
}

# ---- Interaction scatter ----
# Each interaction contrast as its two component responses: the reference line's
# stimulus response on x, the test line's on y.
inter_reg <- contrast_registry |>
  dplyr::filter(type == "interaction")

inter_data <- purrr::pmap_dfr(
  list(inter_reg$name, as.character(inter_reg$experiment), as.character(inter_reg$line),
       as.character(inter_reg$stim), as.character(inter_reg$ref_line),
       inter_reg$label, inter_reg$min_n),
  function(nm, ex, ln, st, ref, lab, mn) {

    test_cn <- paste0(ln, "_", st, "_vs_Mock")
    ctl_cn  <- paste0(ref, "_", st, "_vs_Mock")

    # Component within-line contrasts are looked up by name; skip if either
    # is absent from the registry.
    if (!all(c(nm, test_cn, ctl_cn) %in% names(tt[[ex]]))) {
      message("Interaction ", nm, ": missing component contrast - skipping")
      return(tibble::tibble())
    }

    int  <- tt[[ex]][[nm]]
    test <- tt[[ex]][[test_cn]]
    ctl  <- tt[[ex]][[ctl_cn]]

    tibble::tibble(
      contrast   = nm,
      experiment = ex,
      label      = lab,
      line       = ln,
      stim       = st,
      ref_line   = ref,
      min_n      = mn,
      ensembl_id = int$ensembl_id,
      SYMBOL     = int$SYMBOL,
      x          = ctl$logFC[match(int$ensembl_id, ctl$ensembl_id)],
      y          = test$logFC[match(int$ensembl_id, test$ensembl_id)],
      int_lfc    = int$logFC,
      int_padj   = int$adj.P.Val
    ) |>
      dplyr::filter(!is.na(x), !is.na(y)) |>
      dplyr::mutate(
        sig       = is_de(int_lfc, int_padj),
        direction = dplyr::case_when(
          sig & int_lfc > 0 ~ paste0("Higher in ", ln),
          sig & int_lfc < 0 ~ paste0("Lower in ", ln),
          TRUE              ~ "NS"
        )
      )
  }
)

if (nrow(inter_data) > 0) {

  test_lines <- unique(inter_data$line)
  dir_cols <- c(
    setNames(rep("#D62728", length(test_lines)), paste0("Higher in ", test_lines)),
    setNames(rep("#1F77B4", length(test_lines)), paste0("Lower in ",  test_lines)),
    NS = "grey75"
  )

  sub_i <- "interaction"

  for (nm in unique(inter_data$contrast)) {

    df  <- dplyr::filter(inter_data, contrast == nm)
    lim <- safe_range(c(df$x, df$y))
    r   <- safe_cor(df$x, df$y)

    lab_df <- df |>
      dplyr::filter(sig, !is.na(SYMBOL), SYMBOL != "") |>
      dplyr::slice_min(int_padj, n = 15, with_ties = FALSE)

    sub_txt <- sprintf("%s   |   %d genes differ (adj.P.Val < %s, |log2FC| > %s)",
                       if (is.na(r)) "r = NA" else sprintf("r = %.2f", r),
                       sum(df$sig), P_CUT, LFC_CUT)

    p <- ggplot(df, aes(x, y)) +
      geom_hline(yintercept = 0, colour = "grey88") +
      geom_vline(xintercept = 0, colour = "grey88") +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
      geom_point(aes(colour = direction), size = 1.2, alpha = 0.6) +
      ggrepel::geom_text_repel(
        data = lab_df, aes(label = SYMBOL),
        size = 2.8, max.overlaps = 20, min.segment.length = 0,
        segment.colour = "grey60", segment.size = 0.3) +
      scale_colour_manual(values = dir_cols, name = NULL, drop = FALSE) +
      coord_fixed(xlim = lim, ylim = lim) +
      labs(title = df$label[1], subtitle = sub_txt,
           caption = if (df$min_n[1] == 1) "Unreplicated (n = 1) - exploratory" else NULL,
           x = paste0(df$ref_line[1], ": ", df$stim[1], " vs Mock (log2FC)"),
           y = paste0(df$line[1], ": ", df$stim[1], " vs Mock (log2FC)")) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            plot.caption = element_text(colour = "#B22222", hjust = 0))

    save_fig(p, paste0("interaction_", nm), width = 6.5, height = 6.5, subdir = sub_i)
  }

  # Overview grid: reference line x stimulus, no gene labels
  ov <- inter_data |>
    dplyr::mutate(stim     = factor(stim, levels = STIM_LEVELS),
                  ref_line = factor(ref_line, levels = LINE_LEVELS))
  lim <- safe_range(c(ov$x, ov$y))

  p_all <- ggplot(ov, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
    geom_point(aes(colour = direction), size = 0.9, alpha = 0.55) +
    facet_grid(ref_line ~ stim, drop = FALSE) +
    scale_colour_manual(values = dir_cols, name = NULL, drop = FALSE) +
    coord_fixed(xlim = lim, ylim = lim) +
    labs(title = "Test line stimulus response vs control line response",
         subtitle = "Off-diagonal genes drive the interaction contrasts",
         x = "Control line: stimulus vs Mock (log2FC)",
         y = "Test line: stimulus vs Mock (log2FC)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA))

  save_fig(p_all, "interaction_overview", width = 11, height = 8, subdir = sub_i)
}

# =============================================================================
# DE overlap
# =============================================================================
sub_o <- file.path("de", "overlap")

# ---- Venn: DE genes per stimulus across the spinal lines ----
venn_lines <- c("C9", "Ctrl3", "Ctrl2")
venn_stims <- c("IFNb", "IFNg", "WNV")

wb <- openxlsx::createWorkbook()

for (st in venn_stims) {

  cnames <- setNames(paste0(venn_lines, "_", st, "_vs_Mock"), venn_lines)
  cnames <- cnames[cnames %in% names(tt$spinal)]
  if (length(cnames) < 2) next

  sig_sets <- lapply(cnames, function(cn) de_genes(tt$spinal[[cn]]))

  p <- ggvenn::ggvenn(sig_sets, fill_color = unname(line_cols[names(sig_sets)]),
                      stroke_size = 0.5, set_name_size = 5, text_size = 4) +
    labs(title = paste0(st, " vs Mock: DE genes by line"))

  save_fig(p, paste0("venn_", st), width = 7, height = 6, subdir = sub_o)

  # Per-gene membership plus each line's stats
  all_genes <- sort(unique(unlist(sig_sets, use.names = FALSE)))
  if (length(all_genes) == 0) next

  memb <- vapply(all_genes, function(g) {
    paste(names(sig_sets)[vapply(sig_sets, function(s) g %in% s, logical(1))],
          collapse = ";")
  }, character(1))

  out_df <- tibble::tibble(ensembl_id = all_genes, membership = unname(memb)) |>
    dplyr::left_join(
      tt$spinal[[cnames[[1]]]] |> dplyr::select(ensembl_id, SYMBOL, entrez_id),
      by = "ensembl_id"
    ) |>
    dplyr::relocate(SYMBOL, entrez_id, .after = ensembl_id)

  for (ln in names(cnames)) {
    stats_ln <- tt$spinal[[cnames[[ln]]]] |>
      dplyr::select(ensembl_id, logFC, adj.P.Val) |>
      dplyr::rename(!!paste0(ln, "_logFC") := logFC,
                    !!paste0(ln, "_adj.P.Val") := adj.P.Val)
    out_df <- dplyr::left_join(out_df, stats_ln, by = "ensembl_id")
  }

  openxlsx::addWorksheet(wb, st)
  openxlsx::writeData(wb, st, out_df)
}

if (length(openxlsx::sheets(wb)) > 0) {
  ensure_dir(dir_res)
  openxlsx::saveWorkbook(wb, file.path(dir_res, "DE_genes_venn.xlsx"), overwrite = TRUE)
}

# ---- UpSet: overlap across a whole contrast family ----
# UpSet rather than Venn: ggvenn caps at four sets and the genotype family has
# eight. Up and down are drawn separately.
upset_family <- function(exp_name, contrast_type, direction, file_name) {

  cns <- contrast_registry$name[contrast_registry$type == contrast_type &
                                contrast_registry$experiment == exp_name]
  cns <- intersect(cns, names(tt[[exp_name]]))
  if (length(cns) < 2) return(invisible(NULL))

  sets <- lapply(cns, function(cn) de_genes(tt[[exp_name]][[cn]], direction = direction))
  names(sets) <- reg_of(cns)$label

  # fromList() misbehaves on empty elements
  sets <- sets[lengths(sets) > 0]

  if (length(sets) < 2) {
    save_fig(empty_plot(file_name,
                        note = sprintf("Fewer than two contrasts have %s genes",
                                       tolower(direction))),
             file_name, width = 8, height = 5, subdir = sub_o)
    return(invisible(NULL))
  }

  save_fig(
    function() {
      UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets),
                    order.by = "freq", nintersects = 30,
                    mainbar.y.label = sprintf("Shared %s genes", tolower(direction)),
                    sets.x.label = "DE genes per contrast",
                    text.scale = c(1.3, 1.1, 1.1, 1, 1.1, 1))
    },
    file_name, width = 12, height = 7, subdir = sub_o)
}

upset_family("spinal", "genotype",    "Up",   "upset_genotype_up")
upset_family("spinal", "genotype",    "Down", "upset_genotype_down")
upset_family("spinal", "within_line", "Up",   "upset_within_line_up")
upset_family("spinal", "within_line", "Down", "upset_within_line_down")

# =============================================================================
# Gene expression panels
# =============================================================================
# One faceted figure per panel per experiment.
GENE_PANELS <- list(
  ISG_core      = c("ISG15", "IFIT1", "IFIT3", "MX1", "OAS1", "RSAD2", "IFI6",
                    "STAT1", "IRF7", "BST2", "USP18", "OASL"),
  IFNg_response = c("GBP1", "GBP2", "CXCL9", "CXCL10", "CXCL11", "IDO1", "HLA-DRA",
                    "TAP1", "PSMB9", "IRF1", "SOCS1", "STAT2"),
  Inflammatory  = c("IL6", "IL1B", "TNF", "CCL2", "CCL5", "CXCL8", "NFKB1",
                    "NFKBIA", "TLR3", "TLR4"),
  ALS_C9        = c("C9orf72", "TARDBP", "SOD1", "FUS", "STMN2", "UNC13A", "NEFL",
                    "NEFH", "MAPT", "GFAP", "AIF1", "TREM2"),
  Neural        = c("MAP2", "TUBB3", "SYP", "RBFOX3", "GFAP", "S100B", "OLIG2",
                    "MBP", "SOX2", "NES", "CHAT", "ISL1")
)

# Bar of group means with SEM, individual samples jittered over it.
gene_panel_plot <- function(d, title, subtitle = NULL, ncol = 4) {

  if (nrow(d) == 0) return(empty_plot(title, subtitle, "No panel gene detected"))

  s <- gene_expr_summary(d)

  ggplot(s, aes(condition, mean_expr, fill = condition)) +
    geom_col(width = 0.75, colour = "black", linewidth = 0.25) +
    geom_errorbar(aes(ymin = mean_expr - sem_expr, ymax = mean_expr + sem_expr),
                  width = 0.25, na.rm = TRUE) +
    geom_point(data = d, aes(condition, logcpm),
               position = position_jitter(width = 0.15, height = 0),
               shape = 21, size = 1.8, colour = "black", fill = "white",
               stroke = 0.4, inherit.aes = FALSE) +
    facet_wrap(~ SYMBOL, scales = "free_y", ncol = ncol) +
    scale_fill_manual(values = cond_cols, guide = "none") +
    labs(x = NULL, y = "log2 CPM", title = title, subtitle = subtitle) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

for (e in EXPERIMENTS) {

  lc  <- logcpm[[e]]
  sm  <- sample_meta[[e]]
  map <- symbol_map(tt[[e]])
  sub_p <- file.path("panels", e)

  for (set_name in names(GENE_PANELS)) {

    want <- GENE_PANELS[[set_name]]
    d    <- gene_expr_long(lc, sm, map, want)
    got  <- unique(d$SYMBOL)
    lost <- setdiff(want, got)

    if (length(lost) > 0) {
      message(e, " / ", set_name, ": not detected - ", paste(lost, collapse = ", "))
    }

    fd <- facet_dims(max(1, length(got)), ncol = 4, w_per = 2.6, h_per = 2.4)

    save_fig(
      gene_panel_plot(d, paste0(set_name, " - ", e),
                      sprintf("%d of %d genes detected", length(got), length(want)),
                      ncol = fd$ncol),
      paste0("panel_", set_name), width = fd$width, height = fd$height, subdir = sub_p)
  }

  # Data-driven panel: genes DE in the most contrasts, ties broken on max |logFC|
  hits <- purrr::map_dfr(contrasts_present(e), function(cn) {
    df <- tt[[e]][[cn]]
    keep <- is_de(df$logFC, df$adj.P.Val)
    tibble::tibble(SYMBOL = df$SYMBOL[keep], lfc = abs(df$logFC[keep]))
  }) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::group_by(SYMBOL) |>
    dplyr::summarise(n_contrasts = dplyr::n(), max_lfc = max(lfc), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_contrasts), dplyr::desc(max_lfc)) |>
    dplyr::slice_head(n = 12)

  d  <- gene_expr_long(lc, sm, map, hits$SYMBOL)
  fd <- facet_dims(max(1, dplyr::n_distinct(d$SYMBOL)), ncol = 4, w_per = 2.6, h_per = 2.4)

  save_fig(
    gene_panel_plot(d, paste0("Most recurrent DE genes - ", e),
                    "Ranked by number of contrasts, then max |log2FC|",
                    ncol = fd$ncol),
    "panel_top_recurrent", width = fd$width, height = fd$height, subdir = sub_p)
}

message("5_Graphs.R complete")
