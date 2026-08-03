# ---- Library ----
library(fgsea)
library(ggvenn)
library(ggplot2)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(grid)
library(pheatmap)
library(openxlsx)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("projects/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Data load ----
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

# logcpm is voom's E matrix (log2 counts per million)
tt          <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)), EXPERIMENTS)
logcpm      <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("logcpm_", e), dir = dir_rds)), EXPERIMENTS)
sample_meta <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("sample_meta_", e), dir = dir_rds)), EXPERIMENTS)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)

# Registry row for a contrast name
reg_of <- function(cn) contrast_registry[match(cn, contrast_registry$name), ]

# ---- PCA ----
make_pca <- function(exp_name, n_top = 500) {

  lc <- logcpm[[exp_name]]
  sm <- sample_meta[[exp_name]]

  gene_vars <- apply(lc, 1, var, na.rm = TRUE)
  top <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(n_top, length(gene_vars)))]

  pca <- stats::prcomp(t(lc[top, , drop = FALSE]), scale. = FALSE)
  percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2))

  pca_data <- tibble::tibble(
    sample = rownames(pca$x),
    PC1    = pca$x[, 1],
    PC2    = pca$x[, 2]
  ) |>
    dplyr::left_join(dplyr::select(sm, sample, line, stim, condition), by = "sample")

  centroids <- pca_data |>
    dplyr::group_by(condition, line, stim) |>
    dplyr::summarise(cPC1 = mean(PC1), cPC2 = mean(PC2), .groups = "drop")

  pca_data <- dplyr::left_join(pca_data, centroids, by = c("condition", "line", "stim"))

  p <- ggplot(pca_data, aes(PC1, PC2)) +
    # legs from each sample to its group centroid
    geom_segment(aes(xend = cPC1, yend = cPC2, colour = condition),
                 alpha = 0.5, linewidth = 0.4, show.legend = FALSE) +
    geom_point(aes(fill = condition, shape = stim),
               size = 3.2, colour = "black", stroke = 0.5) +
    geom_point(data = centroids, aes(cPC1, cPC2, fill = condition),
               shape = 21, size = 5, colour = "black", stroke = 1.1,
               inherit.aes = FALSE, show.legend = FALSE) +
    scale_fill_manual(values = cond_cols, drop = TRUE) +
    scale_colour_manual(values = cond_cols, drop = TRUE) +
    scale_shape_manual(values = stim_shapes, drop = TRUE) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(title = paste0("PCA - ", exp_name, " (top ", length(top), " variable genes)")) +
    theme_classic()

  ggsave(file.path(dir_fig, paste0("PCA_", exp_name, ".png")), p,
         width = 8, height = 6, dpi = 300)

  p
}

ensure_dir(dir_fig)
pca_plots <- setNames(lapply(EXPERIMENTS, make_pca), EXPERIMENTS)

# ---- Volcano Plots ----
for (e in EXPERIMENTS) {

  outdir <- ensure_dir(dir_fig, "volcano", e)

  for (coef_name in names(tt[[e]])) {

    volcano_df <- tt[[e]][[coef_name]] |>
      dplyr::mutate(
        neglog10_padj = -log10(adj.P.Val),
        direction = dplyr::case_when(
          adj.P.Val < P_CUT & logFC >  LFC_CUT ~ "Up",
          adj.P.Val < P_CUT & logFC < -LFC_CUT ~ "Down",
          TRUE ~ "NS"
        )
      )

    n_up   <- sum(volcano_df$direction == "Up", na.rm = TRUE)
    n_down <- sum(volcano_df$direction == "Down", na.rm = TRUE)

    x_max <- max(volcano_df$logFC, na.rm = TRUE)
    x_min <- min(volcano_df$logFC, na.rm = TRUE)
    y_max <- max(volcano_df$neglog10_padj[is.finite(volcano_df$neglog10_padj)], na.rm = TRUE)

    p <- ggplot(volcano_df, aes(x = logFC, y = neglog10_padj)) +
      geom_point(aes(color = direction), alpha = 0.6, size = 1) +
      geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed") +
      geom_hline(yintercept = -log10(P_CUT), linetype = "dashed") +
      annotate("text", x = x_max, y = y_max, label = paste0("Up: ", n_up),
               hjust = 1, vjust = 1, color = "red") +
      annotate("text", x = x_min, y = y_max, label = paste0("Down: ", n_down),
               hjust = 0, vjust = 1, color = "blue") +
      scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
      theme_classic() +
      labs(
        title    = reg_of(coef_name)$label,
        subtitle = coef_name,
        x        = "log2 fold change",
        y        = "-log10 adjusted p-value"
      )

    ggsave(file.path(outdir, paste0("volcano_", coef_name, ".png")),
           p, width = 7, height = 6, dpi = 300)
  }
}

# ---- Heatmap: most variable genes ----
for (e in EXPERIMENTS) {

  lc <- logcpm[[e]]
  sm <- sample_meta[[e]]

  sym_map <- tt[[e]][[1]] |>
    dplyr::select(ensembl_id, SYMBOL) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)

  gene_vars <- apply(lc, 1, var, na.rm = TRUE)
  top_ens <- names(sort(gene_vars, decreasing = TRUE))[seq_len(min(50, length(gene_vars)))]

  mat_subset <- lc[top_ens, , drop = FALSE]

  # Label rows with SYMBOL, falling back to Ensembl ID
  new_names <- sym_map$SYMBOL[match(sub("\\..*$", "", top_ens), sym_map$ensembl_id)]
  new_names[is.na(new_names) | new_names == ""] <- top_ens[is.na(new_names) | new_names == ""]
  rownames(mat_subset) <- make.unique(new_names)

  # Z-score per gene
  mat_scaled <- t(scale(t(mat_subset)))

  annotation_col <- data.frame(
    line = sm$line,
    stim = sm$stim,
    row.names = sm$sample
  )

  desired_order <- sm |>
    dplyr::arrange(condition, rep) |>
    dplyr::pull(sample)

  mat_scaled     <- mat_scaled[, desired_order, drop = FALSE]
  annotation_col <- annotation_col[desired_order, , drop = FALSE]

  ann_colors <- list(
    line = line_cols[levels(droplevels(sm$line))],
    stim = setNames(
      grDevices::grey.colors(nlevels(droplevels(sm$stim)), start = 0.2, end = 0.9),
      levels(droplevels(sm$stim))
    )
  )

  pheatmap::pheatmap(
    mat_scaled,
    cluster_rows      = FALSE,
    cluster_cols      = FALSE,
    annotation_col    = annotation_col,
    annotation_colors = ann_colors,
    show_rownames     = TRUE,
    show_colnames     = FALSE,
    main     = paste0("Top 50 variable genes - ", e),
    filename = file.path(dir_fig, paste0("heatmap_top50_variable_", e, ".png")),
    width    = 10,
    height   = 9
  )
}

# ---- log2FC heatmap across contrasts ----
top_n_genes <- 250

for (e in EXPERIMENTS) {

  res_use <- tt[[e]]

  res_use_df <- purrr::imap(res_use, function(df, contrast_name) {
    df |> dplyr::mutate(contrast = contrast_name)
  })

  sig_tbl <- dplyr::bind_rows(res_use_df) |>
    dplyr::filter(!is.na(adj.P.Val), !is.na(logFC)) |>
    dplyr::filter(adj.P.Val < P_CUT, abs(logFC) >= LFC_CUT) |>
    dplyr::select(ensembl_id, logFC, adj.P.Val, contrast)

  if (nrow(sig_tbl) == 0) {
    message("No significant genes for ", e, " - skipping log2FC heatmap")
    next
  }

  top_genes <- sig_tbl |>
    dplyr::group_by(ensembl_id) |>
    dplyr::summarise(max_abs_lfc = max(abs(logFC)), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(max_abs_lfc)) |>
    dplyr::slice_head(n = top_n_genes) |>
    dplyr::pull(ensembl_id)

  # Non-significant cells become NA and are drawn black
  lfc_mat <- sapply(res_use_df, function(df) {
    idx <- match(top_genes, df$ensembl_id)
    lfc <- df$logFC[idx]
    p   <- df$adj.P.Val[idx]

    is_sig <- !is.na(p) & !is.na(lfc) & (p < P_CUT) & (abs(lfc) >= LFC_CUT)
    lfc[!is_sig] <- NA_real_
    lfc
  })
  colnames(lfc_mat) <- names(res_use_df)

  reg_e <- contrast_registry |>
    dplyr::filter(experiment == e) |>
    dplyr::arrange(type, line, ref_line, stim)

  lfc_mat <- lfc_mat[, reg_e$name, drop = FALSE]

  annotation_col <- data.frame(type = reg_e$type, row.names = reg_e$label)
  colnames(lfc_mat) <- reg_e$label

  sym_map <- res_use[[1]] |>
    dplyr::select(ensembl_id, SYMBOL) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)

  gene_labels <- sym_map$SYMBOL[match(top_genes, sym_map$ensembl_id)]
  gene_labels[is.na(gene_labels) | gene_labels == ""] <- top_genes[is.na(gene_labels) | gene_labels == ""]
  rownames(lfc_mat) <- make.unique(gene_labels)

  # NA -> 0 for the clustering distance only
  lfc_for_cluster <- lfc_mat
  lfc_for_cluster[is.na(lfc_for_cluster)] <- 0
  row_hc <- stats::hclust(stats::dist(lfc_for_cluster))

  pheatmap::pheatmap(
    lfc_mat,
    cluster_rows      = row_hc,
    cluster_cols      = FALSE,
    annotation_col    = annotation_col,
    annotation_colors = list(type = type_cols[levels(droplevels(reg_e$type))]),
    na_col            = "black",
    show_colnames     = TRUE,
    show_rownames     = FALSE,
    angle_col         = 45,
    main     = paste0("Top ", nrow(lfc_mat), " DE genes, log2 fold change - ", e),
    filename = file.path(dir_fig, paste0("heatmap_log2FC_", e, ".png")),
    width    = 11,
    height   = 10
  )
}

# ---- Bar graphs for a hallmark gene set ----
# EDIT: pick a hallmark set from script 2's GSEA results
hallmark_name <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
symbols <- unique(hallmark_sets[[hallmark_name]])

for (e in EXPERIMENTS) {

  res_use <- tt[[e]]
  lc <- logcpm[[e]]
  sm <- sample_meta[[e]]

  # Set members significant in at least one contrast of this experiment
  sig_ens <- dplyr::bind_rows(res_use) |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL %in% symbols) |>
    dplyr::filter(!is.na(adj.P.Val), !is.na(logFC)) |>
    dplyr::filter(adj.P.Val < P_CUT, abs(logFC) >= LFC_CUT) |>
    dplyr::pull(ensembl_id) |>
    unique()

  if (length(sig_ens) == 0) {
    message("No significant ", hallmark_name, " genes for ", e, " - skipping bar graphs")
    next
  }

  sym_map <- res_use[[1]] |>
    dplyr::select(ensembl_id, SYMBOL) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)

  keep_rows <- sub("\\..*$", "", rownames(lc)) %in% sig_ens
  lc_sub <- lc[keep_rows, , drop = FALSE]

  # Wide to long: one row per gene-sample
  plot_df <- tibble::tibble(
    ensembl_id = rep(sub("\\..*$", "", rownames(lc_sub)), times = ncol(lc_sub)),
    sample     = rep(colnames(lc_sub), each = nrow(lc_sub)),
    logcpm     = as.vector(lc_sub)
  ) |>
    dplyr::left_join(dplyr::select(sm, sample, line, stim, condition), by = "sample") |>
    dplyr::left_join(sym_map, by = "ensembl_id") |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL %in% symbols)

  summary_df <- plot_df |>
    dplyr::group_by(SYMBOL, condition, line, stim) |>
    dplyr::summarise(
      n         = sum(!is.na(logcpm)),
      mean_expr = mean(logcpm, na.rm = TRUE),
      sd_expr   = stats::sd(logcpm, na.rm = TRUE),
      sem_expr  = sd_expr / sqrt(n),
      .groups   = "drop"
    )

  # Bar of group means with SEM, individual samples jittered over it
  make_gene_plot <- function(gene_symbol) {
    df_pts <- dplyr::filter(plot_df, SYMBOL == gene_symbol)
    df_sum <- dplyr::filter(summary_df, SYMBOL == gene_symbol)

    ggplot(df_sum, aes(x = condition, y = mean_expr, fill = condition)) +
      geom_col(width = 0.75, colour = "black") +
      geom_errorbar(aes(ymin = mean_expr - sem_expr, ymax = mean_expr + sem_expr),
                    width = 0.25, na.rm = TRUE) +
      geom_point(
        data = df_pts,
        mapping = aes(x = condition, y = logcpm),
        position = position_jitter(width = 0.15, height = 0),
        shape = 21, size = 2.5, colour = "black", fill = "white", stroke = 0.6
      ) +
      scale_fill_manual(values = cond_cols, guide = "none") +
      labs(x = NULL, y = "log2 CPM", title = gene_symbol) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  }

  sig_symbols <- summary_df |>
    dplyr::filter(!is.na(SYMBOL), SYMBOL != "") |>
    dplyr::distinct(SYMBOL) |>
    dplyr::pull(SYMBOL)

  outdir <- ensure_dir(dir_fig, "bargraphs", e, hallmark_name)

  for (g in sig_symbols) {
    ggsave(file.path(outdir, paste0(g, ".png")), make_gene_plot(g),
           width = 7, height = 4, dpi = 300)
  }

  message(e, ": wrote ", length(sig_symbols), " bar graphs to ", outdir)
}

# ---- Venn diagrams ----
# One 3-set Venn per stimulus, comparing DE genes across the spinal lines
venn_lines <- c("C9", "Ctrl3", "Ctrl2")
venn_stims <- c("IFNb", "IFNg", "WNV")

sig_genes_for <- function(df) {
  df |>
    dplyr::filter(!is.na(adj.P.Val), !is.na(logFC)) |>
    dplyr::filter(adj.P.Val < P_CUT, abs(logFC) >= LFC_CUT) |>
    dplyr::pull(ensembl_id) |>
    unique()
}

outdir <- ensure_dir(dir_fig, "venn")
wb <- openxlsx::createWorkbook()

for (st in venn_stims) {

  cnames <- setNames(paste0(venn_lines, "_", st, "_vs_Mock"), venn_lines)
  cnames <- cnames[cnames %in% names(tt$spinal)]
  if (length(cnames) < 2) next

  sig_sets <- lapply(cnames, function(cn) sig_genes_for(tt$spinal[[cn]]))

  p <- ggvenn::ggvenn(sig_sets, fill_color = unname(line_cols[names(sig_sets)]),
                      stroke_size = 0.5, set_name_size = 5, text_size = 4) +
    labs(title = paste0(st, " vs Mock: DE genes by line"))

  ggsave(file.path(outdir, paste0("venn_", st, ".png")), p,
         width = 7, height = 6, dpi = 300)

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
