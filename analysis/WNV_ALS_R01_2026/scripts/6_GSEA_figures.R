# GSEA figures for all 25 contrasts, from the checkpoints 4_GSEA.R wrote.
# Run from the repo root, after 5_Graphs.R.
#
# fgsea() is never called here. Its p-values are permutation-based and 4_GSEA.R
# sets no seed, so a re-run would not reproduce results/gsea_hallmark_all.tsv.
# Only fgsea::plotEnrichment() is used, which is deterministic given the ranks.

# ---- Library ----
library(fgsea)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(stringr)
library(grid)
library(pheatmap)
library(patchwork)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Parameters ----
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
GMT_PATH <- "genesets/h.all.v2025.1.Hs.symbols.gmt"

# A pathway is drawn on the bubble plot if it is significant in at least this
# many contrasts of the experiment.
BUBBLE_MIN_SIG <- 1L

# Pathways per contrast on the lollipop, and enrichment curves per contrast.
GSEA_TOP_N  <- 15
CURVE_TOP_N <- 5

# Leading-edge heatmaps are built for this fixed list only
LE_SETS <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
)
LE_MAX_GENES <- 60

sig_label <- paste0("padj < ", P_CUT)

# ---- Data load ----
hallmark_sets <- fgsea::gmtPathways(GMT_PATH)

tt          <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)), EXPERIMENTS)
logcpm      <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("logcpm_", e), dir = dir_rds)), EXPERIMENTS)
sample_meta <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("sample_meta_", e), dir = dir_rds)), EXPERIMENTS)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)
gsea_hallmark     <- load_checkpoint("gsea_hallmark", dir = dir_rds)
gsea_als          <- load_checkpoint("gsea_als",      dir = dir_rds)

contrasts_present <- function(e) intersect(contrasts_for(e), names(tt[[e]]))

# ---- Ranked statistics ----
# Rebuilt with the recipe run_gsea() used in 4_GSEA.R so the running-ES curves
# match the checkpointed NES. The block below asserts the two agree.
ranks <- purrr::map(tt, function(tts) {
  purrr::map(tts, rank_stats, gene_col = "SYMBOL", stat_col = "t",
             dedupe_by = "AveExpr", quiet = TRUE)
})

local({
  bad <- character(0)
  for (e in EXPERIMENTS) for (cn in contrasts_present(e)) {
    x <- as.data.frame(tt[[e]][[cn]])
    x <- x[!is.na(x$SYMBOL) & x$SYMBOL != "", ]
    x <- x[order(-x$AveExpr), ]
    x <- x[!duplicated(x$SYMBOL), ]
    if (!identical(sort(setNames(x$t, x$SYMBOL), decreasing = TRUE), ranks[[e]][[cn]]))
      bad <- c(bad, cn)
  }
  if (length(bad) > 0)
    stop("rank_stats() diverges from run_gsea()'s recipe for: ", paste(bad, collapse = ", "))
  message("rank_stats matches run_gsea() for all ",
          sum(lengths(lapply(EXPERIMENTS, contrasts_present))), " contrasts")
})

sub_hall <- file.path("gsea", "hallmark")
sub_als  <- file.path("gsea", "als")

# =============================================================================
# 1. Hallmark bubble: pathway x contrast
# =============================================================================
# Non-significant cells are drawn as faint crosses rather than filtered out, so
# a contrast with nothing significant keeps its column.
make_gsea_bubble <- function(exp_name) {

  lv <- contrast_labels(exp_name)

  df <- gsea_hallmark |>
    dplyr::filter(experiment == exp_name, !is.na(NES), !is.na(pathway)) |>
    dplyr::mutate(sig = !is.na(padj) & padj <= P_CUT) |>
    dplyr::group_by(pathway) |>
    dplyr::filter(sum(sig) >= BUBBLE_MIN_SIG) |>
    dplyr::ungroup()

  if (nrow(df) == 0) {
    return(empty_plot(paste0("Hallmark GSEA - ", exp_name),
                      note = "No pathway significant in any contrast"))
  }

  # Scaffold every retained pathway against every registry contrast so a
  # contrast absent from gsea_hallmark still gets a labelled column.
  # scale_x_discrete(drop = FALSE) does not work here: with
  # facet_grid(scales = "free_x") every facet would show every contrast.
  reg_e <- contrast_registry |>
    dplyr::filter(experiment == exp_name) |>
    dplyr::select(label, type, min_n)

  df <- tidyr::expand_grid(pathway = unique(df$pathway), reg_e) |>
    dplyr::left_join(dplyr::select(df, pathway, label, NES, padj, sig),
                     by = c("pathway", "label")) |>
    dplyr::mutate(sig = !is.na(sig) & sig)

  df <- df |>
    dplyr::mutate(
      label   = factor(label, levels = lv),
      pathway = stats::reorder(sub("^HALLMARK_", "", pathway), NES,
                               FUN = function(z) mean(z, na.rm = TRUE)),
      # Unreplicated groups drawn as triangles; the factor fixes the legend
      # across experiments
      replication = factor(ifelse(min_n == 1, "n = 1", "replicated"),
                           levels = c("replicated", "n = 1"))
    )

  ggplot(df, aes(label, pathway)) +
    geom_point(data = ~ dplyr::filter(.x, !sig), shape = 4, size = 1,
               colour = "grey82", show.legend = FALSE) +
    geom_point(data = ~ dplyr::filter(.x, sig),
               aes(size = -log10(padj), fill = NES, shape = replication),
               colour = "black") +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    scale_x_discrete(position = "top") +
    scale_shape_manual(values = c("n = 1" = 24, "replicated" = 21), drop = FALSE) +
    scale_size(name = "-log10(padj)", range = c(3, 9)) +
    # Symmetric limits put NES 0 at the palette midpoint
    scale_fill_distiller(palette = "Spectral",
                         limits = c(-1, 1) * safe_absmax(df$NES[df$sig])) +
    labs(title = paste0("Hallmark GSEA - ", exp_name),
         subtitle = paste0("grey cross = tested, not significant at ", sig_label),
         x = NULL, y = "Gene Set") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(colour = "black"),
          axis.text.x = element_text(colour = "black", angle = 45, hjust = 0),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

# =============================================================================
# 2. NES heatmap: every pathway x contrast cell drawn
# =============================================================================
make_nes_heatmap <- function(exp_name) {

  lv    <- contrast_labels(exp_name)
  reg_e <- contrast_registry |>
    dplyr::filter(experiment == exp_name) |>
    dplyr::arrange(type, line, ref_line, stim)

  d <- gsea_hallmark |> dplyr::filter(experiment == exp_name)
  if (nrow(d) == 0) return(NULL)

  # complete() guarantees a cell for every pathway x contrast, including sets
  # fgsea dropped on minSize.
  w <- d |>
    dplyr::mutate(label = factor(label, levels = lv)) |>
    dplyr::select(pathway, label, NES, padj) |>
    tidyr::complete(pathway, label)

  to_mat <- function(col) {
    m <- w |>
      dplyr::select(pathway, label, dplyr::all_of(col)) |>
      tidyr::pivot_wider(names_from = label, values_from = dplyr::all_of(col)) |>
      tibble::column_to_rownames("pathway") |>
      as.matrix()
    m[, lv, drop = FALSE]
  }

  m <- to_mat("NES")
  stopifnot(is.numeric(m))   # a duplicate key would give a list-column here

  stars <- w |>
    dplyr::mutate(s = dplyr::case_when(
      is.na(padj)    ~ "",
      padj <= 0.001  ~ "***",
      padj <= 0.01   ~ "**",
      padj <= P_CUT  ~ "*",
      TRUE           ~ "")) |>
    dplyr::select(pathway, label, s) |>
    tidyr::pivot_wider(names_from = label, values_from = s) |>
    tibble::column_to_rownames("pathway") |>
    as.matrix()
  stars <- stars[rownames(m), lv, drop = FALSE]

  rownames(m) <- rownames(stars) <- abbrev_pathway(rownames(m), 34)

  amax <- safe_absmax(m)

  # NA -> 0 for the clustering distance only
  mc <- m
  mc[is.na(mc)] <- 0

  pheatmap::pheatmap(
    m, silent = TRUE,
    cluster_rows = if (heatmap_ready(mc)) stats::hclust(stats::dist(mc)) else FALSE,
    cluster_cols = FALSE,
    display_numbers = stars, fontsize_number = 7, number_color = "black",
    color  = grDevices::colorRampPalette(
               rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    breaks = seq(-amax, amax, length.out = 101),
    na_col = "grey90",
    annotation_col = data.frame(
      type  = reg_e$type,
      min_n = ifelse(reg_e$min_n == 1, "n = 1", "replicated"),
      row.names = reg_e$label
    ),
    annotation_colors = list(
      type  = type_cols[levels(droplevels(reg_e$type))],
      min_n = c("n = 1" = "#B22222", "replicated" = "grey85")
    ),
    angle_col = 45,
    main = sprintf("Hallmark NES - %s   (* padj<%s, ** <0.01, *** <0.001)",
                   exp_name, P_CUT)
  )
}

# =============================================================================
# 3. Per-contrast lollipop
# =============================================================================
make_gsea_lollipop <- function(exp_name, contrast_name) {

  r <- reg_of(contrast_name)
  d <- gsea_hallmark |>
    dplyr::filter(experiment == exp_name, contrast == contrast_name, !is.na(NES))

  if (nrow(d) == 0) return(empty_plot(r$label, contrast_name, "No GSEA result"))

  sig      <- dplyr::filter(d, !is.na(padj), padj <= P_CUT)
  fallback <- nrow(sig) == 0

  top <- if (fallback) {
    dplyr::slice_min(dplyr::filter(d, !is.na(pval)), pval, n = GSEA_TOP_N, with_ties = FALSE)
  } else {
    dplyr::slice_max(sig, abs(NES), n = GSEA_TOP_N, with_ties = FALSE)
  }

  if (nrow(top) == 0) return(empty_plot(r$label, contrast_name, "No GSEA result"))

  top <- top |>
    dplyr::arrange(NES) |>
    dplyr::mutate(clean     = abbrev_pathway(pathway, 30),
                  direction = ifelse(NES > 0, "Enriched", "Suppressed"))
  top$clean <- factor(top$clean, levels = unique(top$clean))

  # Arrows sit in the space added below the first row; symmetric x limits put
  # one label in each half panel.
  xm <- safe_absmax(top$NES) * 1.15

  up_lab <- if (!is.na(r$ref_line)) paste0("Higher in ", r$line) else
            paste0("Higher in ", r$stim)
  dn_lab <- if (!is.na(r$ref_line)) paste0("Higher in ", r$ref_line) else "Higher in Mock"

  sub_txt <- sprintf("%d of %d sets at %s%s", nrow(sig), nrow(d), sig_label,
                     if (isTRUE(r$min_n == 1)) "   |   unreplicated (n = 1)" else "")

  ggplot(top, aes(NES, clean)) +
    geom_segment(aes(x = 0, xend = NES, yend = clean),
                 colour = "grey60", linewidth = 0.4) +
    geom_point(aes(size = -log10(padj), fill = direction), shape = 21, stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    annotate("segment", x = 0, xend = xm, y = -0.2, yend = -0.2,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             colour = "grey30", linewidth = 0.9) +
    annotate("segment", x = 0, xend = -xm, y = -0.2, yend = -0.2,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed"),
             colour = "grey30", linewidth = 0.9) +
    annotate("text", x = xm, y = -1.0, label = up_lab, hjust = 1, size = 3.2) +
    annotate("text", x = -xm, y = -1.0, label = dn_lab, hjust = 0, size = 3.2) +
    scale_fill_manual(values = c(Enriched = "#D73027", Suppressed = "#4575B4"),
                      guide = "none") +
    scale_size_continuous(name = "-log10(padj)", range = c(3, 8)) +
    scale_x_continuous(limits = c(-xm, xm), expand = expansion(mult = 0.03)) +
    scale_y_discrete(expand = expansion(add = c(3.0, 0.6))) +
    labs(title = r$label, subtitle = sub_txt,
         caption = if (fallback)
           sprintf("Nothing reaches %s; showing top %d by nominal p", sig_label, GSEA_TOP_N)
           else NULL,
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          text = element_text(colour = "black"),
          plot.title = element_text(face = "bold"),
          plot.caption = element_text(colour = "#B22222", hjust = 0))
}

# =============================================================================
# 4. Running-ES curves, one patchwork per contrast
# =============================================================================
make_curve_panel <- function(exp_name, contrast_name) {

  r  <- reg_of(contrast_name)
  rk <- ranks[[exp_name]][[contrast_name]]
  d  <- gsea_hallmark |>
    dplyr::filter(experiment == exp_name, contrast == contrast_name, !is.na(padj))

  top <- dplyr::slice_min(d, padj, n = CURVE_TOP_N, with_ties = FALSE)

  if (length(rk) == 0 || nrow(top) == 0) {
    return(empty_plot(r$label, contrast_name, "No ranked statistics / no GSEA result"))
  }

  ps <- purrr::pmap(list(top$pathway, top$NES, top$padj), function(pw, nes, padj) {
    genes <- intersect(hallmark_sets[[pw]], names(rk))
    # plotEnrichment errors on an empty intersection
    if (length(genes) < 2) return(empty_plot(abbrev_pathway(pw, 30)))
    fgsea::plotEnrichment(genes, rk) +
      labs(title = abbrev_pathway(pw, 30),
           subtitle = sprintf("NES %.2f   padj %.2g", nes, padj),
           x = "Rank in ordered gene list", y = "Enrichment score") +
      theme(plot.title    = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 8))
  })

  patchwork::wrap_plots(ps, ncol = min(3, length(ps))) +
    patchwork::plot_annotation(
      title = r$label,
      subtitle = sprintf("%s   |   top %d sets by padj%s", contrast_name, nrow(top),
                         if (isTRUE(r$min_n == 1)) "   |   unreplicated (n = 1)" else ""))
}

# =============================================================================
# 5. Leading-edge heatmaps
# =============================================================================
# Leading-edge genes x samples, z-scored. The genes x contrasts view would
# duplicate the NES heatmap and overview/heatmap_log2FC_<e>.
make_le_heatmap <- function(exp_name, pw) {

  le <- gsea_hallmark |>
    dplyr::filter(experiment == exp_name, pathway == pw,
                  !is.na(padj), padj <= P_CUT) |>
    dplyr::pull(leadingEdge) |>
    split_leading_edge()

  if (length(le) < 2) {
    message(exp_name, " / ", pw, ": leading edge has <2 genes - skipping")
    return(NULL)
  }

  lc  <- logcpm[[exp_name]]
  sm  <- sample_meta[[exp_name]]
  map <- symbol_map(tt[[exp_name]])

  d <- gene_expr_long(lc, sm, map, le)
  if (nrow(d) == 0) return(NULL)

  m <- d |>
    dplyr::select(SYMBOL, sample, logcpm) |>
    tidyr::pivot_wider(names_from = sample, values_from = logcpm) |>
    tibble::column_to_rownames("SYMBOL") |>
    as.matrix()

  # Drop constant rows (all-NaN after scaling), keep the most variable members,
  # then z-score.
  m <- drop_constant_rows(m)
  if (nrow(m) > LE_MAX_GENES) {
    v <- matrixStats::rowVars(m)
    m <- m[utils::head(order(v, decreasing = TRUE), LE_MAX_GENES), , drop = FALSE]
  }
  m <- zscore_rows(m)

  if (!heatmap_ready(m)) {
    message(exp_name, " / ", pw, ": too few usable genes - skipping")
    return(NULL)
  }

  ord <- sm |> dplyr::arrange(condition, rep) |> dplyr::pull(sample)
  ord <- intersect(ord, colnames(m))
  m   <- m[, ord, drop = FALSE]

  ann <- data.frame(line = sm$line, stim = sm$stim, row.names = sm$sample)[ord, , drop = FALSE]

  n_sig <- gsea_hallmark |>
    dplyr::filter(experiment == exp_name, pathway == pw, !is.na(padj), padj <= P_CUT) |>
    nrow()

  pheatmap::pheatmap(
    m, silent = TRUE,
    cluster_rows = TRUE, cluster_cols = FALSE,
    annotation_col = ann,
    annotation_colors = list(
      line = line_cols[levels(droplevels(sm$line))],
      stim = setNames(grDevices::grey.colors(nlevels(droplevels(sm$stim)),
                                             start = 0.2, end = 0.9),
                      levels(droplevels(sm$stim)))
    ),
    show_colnames = FALSE,
    main = sprintf("%s\nleading edge, %d genes, union over %d significant contrast(s) - %s",
                   sub("^HALLMARK_", "", pw), nrow(m), n_sig, exp_name)
  )
}

# =============================================================================
# 6. Custom (ALS) gene set NES
# =============================================================================
make_als_bar <- function(exp_name) {

  lv <- contrast_labels(exp_name)

  # Scaffolded off the registry so a contrast with no row in gsea_als still gets
  # a labelled slot, as in make_gsea_bubble().
  df <- contrast_registry |>
    dplyr::filter(experiment == exp_name) |>
    dplyr::select(contrast = name, label, type, min_n) |>
    dplyr::left_join(
      dplyr::select(gsea_als, contrast, term, NES, p.adjust),
      by = "contrast") |>
    dplyr::mutate(
      label       = factor(label, levels = lv),
      significant = dplyr::case_when(
        is.na(p.adjust)      ~ "no result",
        p.adjust <= P_CUT    ~ sig_label,
        TRUE                 ~ "ns"),
      significant = factor(significant, levels = c(sig_label, "ns", "no result")),
      NES_plot    = ifelse(is.na(NES), 0, NES)
    )

  if (nrow(df) == 0) {
    return(empty_plot(paste0("Custom gene set - ", exp_name),
                      note = "No result for this experiment"))
  }

  term_name <- unique(df$term[!is.na(df$term)])
  title_txt <- if (length(term_name) == 1) {
    str_to_sentence(str_replace_all(term_name, "_", " "))
  } else {
    "Custom gene set"
  }

  ggplot(df, aes(label, NES_plot, fill = significant)) +
    geom_col(colour = "black", width = 0.7) +
    # Unreplicated contrasts flagged at the free end of the bar
    geom_text(aes(label = ifelse(min_n == 1, "n = 1", ""),
                  vjust = ifelse(NES_plot >= 0, -0.4, 1.3)),
              size = 2.6, colour = "#B22222") +
    geom_hline(yintercept = 0) +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = setNames(c("#D62728", "grey80", "white"),
                                        c(sig_label, "ns", "no result")),
                      drop = FALSE, name = NULL) +
    labs(x = NULL, y = "NES", title = title_txt, subtitle = exp_name) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
          strip.background = element_rect(fill = "grey92", colour = NA))
}

# =============================================================================
# Draw
# =============================================================================
for (e in EXPERIMENTS) {

  cns <- contrasts_present(e)
  n_x <- length(contrast_labels(e))
  w   <- max(7, 2.2 + n_x * 0.9)

  # --- overview ---
  save_fig(make_gsea_bubble(e), paste0("bubble_", e),
           width = w, height = 9, subdir = sub_hall)

  ph <- make_nes_heatmap(e)
  if (!is.null(ph)) {
    save_fig(ph, paste0("nes_heatmap_", e),
             width = max(8, 4.5 + n_x * 0.55), height = 11, subdir = sub_hall)
  }

  save_fig(make_als_bar(e), paste0("als_NES_", e),
           width = w, height = 5, subdir = sub_als)

  # --- per contrast ---
  for (cn in cns) {
    save_fig(make_gsea_lollipop(e, cn), paste0("top_pathways_", cn),
             width = 8, height = 7, subdir = file.path(sub_hall, "lollipop", e))

    save_fig(make_curve_panel(e, cn), paste0("curves_", cn),
             width = 12, height = 7, subdir = file.path(sub_hall, "curves", e))
  }

  # --- leading edge ---
  for (pw in LE_SETS) {
    ph <- make_le_heatmap(e, pw)
    if (!is.null(ph)) {
      save_fig(ph, paste0("le_", sub("^HALLMARK_", "", pw)),
               width = max(8, 2 + ncol(logcpm[[e]]) * 0.22), height = 11,
               subdir = file.path(sub_hall, "leading_edge", e))
    }
  }
}

message("6_GSEA_figures.R complete")
