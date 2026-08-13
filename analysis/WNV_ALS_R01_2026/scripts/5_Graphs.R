# ---- Library ----
library(fgsea)
library(ggvenn)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(grid)
library(pheatmap)
library(openxlsx)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# =============================================================================
# Original Figures
# =============================================================================
# Objects in data/fig4_r2sdhf/ come from analysis/R2SDHF/data/r_objects
#+ FIGURE 4 - R2SDHF WNV infection timecourse
fig4_dir <- file.path(PROJ, "data", "fig4_r2sdhf")

FIG4_NES_CUT  <- 1.5
FIG4_PVAL_CUT <- 0.05
FIG4_N_LABELS <- 25

# Fonts, pt
FIG4_AXIS_TITLE <- 20
FIG4_AXIS_TEXT  <- 16
FIG4_TITLE_SIZE <- 22
# Direction arrow labels, pt
FIG4_ARROW_TEXT <- 14
# Gene labels on the scatter, pt
FIG4_GENE_TEXT  <- 14

# Pathway labels truncated past this many characters
FIG4_LABEL_MAX <- 30

# Panel height per pathway row, relative to panel width
FIG4_ROW_ASPECT <- 0.085
# Bubble size in mm; -log10(padj) limits shared by all three panels
FIG4_SIZE_RANGE  <- c(3, 9)
FIG4_SIZE_LIMITS <- c(0, 50)
FIG4_SIZE_BREAKS <- c(1, 25, 50)

# Scatter gene list
FIG4_INNATE_SETS <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT"
)

# Bubble plot pathways
FIG4_BUBBLE_SETS <- c(
  FIG4_INNATE_SETS,
  "HALLMARK_P53_PATHWAY",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)

fig4_times <- c(TX12 = "12", TX24 = "24", TX48 = "48")

fig4_gsea <- setNames(lapply(names(fig4_times), function(tp) {
  readRDS(file.path(fig4_dir, sprintf("fgsea_%s_vs_Mock(R2SDHF).rds", tp)))
}), names(fig4_times))

fig4_res <- setNames(lapply(names(fig4_times), function(tp) {
  readRDS(file.path(fig4_dir, sprintf("res_%s(R2SDHF).rds", tp)))
}), names(fig4_times))

# ---- Fig 4A-C: hallmark GSEA, WNV vs Mock
fig4_bubble_data <- function(gsea_res) {
  gsea_res |>
    dplyr::filter(pathway %in% FIG4_BUBBLE_SETS,
                  abs(NES) > FIG4_NES_CUT, pval < FIG4_PVAL_CUT) |>
    dplyr::mutate(
      pathway_clean  = abbrev_pathway(pathway, FIG4_LABEL_MAX),
      direction      = ifelse(NES > 0, "Enriched", "Suppressed"),
      neg_log10_padj = -log10(padj)
    ) |>
    dplyr::arrange(NES)
}

fig4_bubble <- function(d, time_label, size_scale) {

  if (nrow(d) == 0) return(NULL)
  d$pathway_clean <- factor(d$pathway_clean, levels = d$pathway_clean)

  # Arrows sit in the space added below the first row; symmetric x limits put
  # one label in each half panel
  xm    <- max(abs(d$NES)) * 1.15
  y_arr <- -0.2
  y_lab <- -1.0

  ggplot(d, aes(NES, pathway_clean)) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_clean),
                 colour = "grey60", linewidth = 0.4) +
    geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    annotate("segment", x = 0, xend = xm, y = y_arr, yend = y_arr,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             colour = "grey30", linewidth = 1.1) +
    annotate("segment", x = 0, xend = -xm, y = y_arr, yend = y_arr,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             colour = "grey30", linewidth = 1.1) +
    annotate("text", x = xm, y = y_lab, label = "Higher in WNV",
             hjust = 1, colour = "black", size = FIG4_ARROW_TEXT / .pt) +
    annotate("text", x = -xm, y = y_lab, label = "Higher in Mock",
             hjust = 0, colour = "black", size = FIG4_ARROW_TEXT / .pt) +
    scale_fill_manual(values = c(Enriched = "#D73027", Suppressed = "#4575B4"),
                      guide = "none") +
    scale_size_continuous(range  = FIG4_SIZE_RANGE,
                          limits = size_scale$limits,
                          breaks = size_scale$breaks,
                          name   = "-log10(p-value)") +
    scale_x_continuous(limits = c(-xm, xm), expand = expansion(mult = 0.03)) +
    scale_y_discrete(expand = expansion(add = c(3.0, 0.6))) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = sprintf("WNV vs. Mock (%s hpi)", time_label)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          aspect.ratio = FIG4_ROW_ASPECT * (nrow(d) + 3),
          text         = element_text(colour = "black"),
          plot.title   = element_text(face = "bold", size = FIG4_TITLE_SIZE, colour = "black"),
          axis.title   = element_text(size = FIG4_AXIS_TITLE, colour = "black"),
          axis.text    = element_text(size = FIG4_AXIS_TEXT, colour = "black"),
          legend.text  = element_text(size = FIG4_AXIS_TEXT, colour = "black"),
          legend.title = element_text(size = FIG4_AXIS_TEXT, colour = "black"))
}

fig4_bubble_dfs <- lapply(fig4_gsea, fig4_bubble_data)

fig4_size <- list(limits = FIG4_SIZE_LIMITS, breaks = FIG4_SIZE_BREAKS)

fig4_bubbles <- purrr::imap(fig4_bubble_dfs, function(d, tp) {
  fig4_bubble(d, fig4_times[[tp]], fig4_size)
})

fig4_bubbles$TX12
fig4_bubbles$TX24
fig4_bubbles$TX48

# ---- Fig 4D: innate immunity genes at 48 h
# Leading-edge genes of the innate sets, 48 h GSEA
fig4_targets <- fig4_gsea$TX48 |>
  dplyr::filter(pathway %in% FIG4_INNATE_SETS) |>
  dplyr::pull(leadingEdge) |>
  unlist() |>
  unique()

fig4_volc <- fig4_res$TX48 |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  dplyr::mutate(
    FC        = 2^logFC,
    log10_p   = log10(P.Value),
    is_target = gene %in% fig4_targets,
    direction = dplyr::case_when(
      is_target & logFC > 0 ~ "Up",
      is_target & logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  )

# y limit set by the least significant highlighted gene
fig4_yfloor <- min(fig4_volc$log10_p[fig4_volc$is_target], na.rm = TRUE) * 1.2

# Band above the axis holding the direction arrows
fig4_ytop <- abs(fig4_yfloor) * 0.18

fig4_scatter <- ggplot(fig4_volc, aes(FC, log10_p)) +
  geom_point(data = function(d) dplyr::filter(d, !is_target),
             colour = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = function(d) dplyr::filter(d, is_target),
             aes(colour = direction), size = 1.8, alpha = 0.9) +
  # geom text size is mm, /.pt converts from pt
  geom_text_repel(data = function(d) dplyr::slice_min(dplyr::filter(d, is_target),
                                                      P.Value, n = FIG4_N_LABELS),
                  aes(label = gene), size = FIG4_GENE_TEXT / .pt, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.4,
                  ylim = c(fig4_yfloor, 0)) +
  geom_vline(xintercept = 1, linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  annotate("segment", x = 1, xend = 40, y = fig4_ytop * 0.3, yend = fig4_ytop * 0.3,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           colour = "grey30", linewidth = 1.1) +
  annotate("segment", x = 1, xend = 1/40, y = fig4_ytop * 0.3, yend = fig4_ytop * 0.3,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           colour = "grey30", linewidth = 1.1) +
  annotate("text", x = 40, y = fig4_ytop * 0.78, label = "Higher in WNV",
           hjust = 1, colour = "black", size = FIG4_ARROW_TEXT / .pt) +
  annotate("text", x = 1/40, y = fig4_ytop * 0.78, label = "Higher in Mock",
           hjust = 0, colour = "black", size = FIG4_ARROW_TEXT / .pt) +
  scale_colour_manual(values = c(Up = "red", Down = "black"), guide = "none") +
  scale_x_continuous(trans = "log2", breaks = 2^(-5:5),
                     labels = c("-32", "-16", "-8", "-4", "-2", "1",
                                "+2", "+4", "+8", "+16", "+32")) +
  coord_cartesian(xlim = c(1/45, 45), ylim = c(fig4_yfloor, fig4_ytop)) +
  labs(x = "FC: WNV vs. Mock", y = "log10(p-value)",
       title = "Innate Immunity Genes (48 hpi)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        text       = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = FIG4_TITLE_SIZE, colour = "black"),
        axis.title = element_text(size = FIG4_AXIS_TITLE, colour = "black"),
        axis.text  = element_text(size = FIG4_AXIS_TEXT, colour = "black"))

fig4_scatter

# ---- Save Figure 4 -
purrr::iwalk(fig4_bubbles, function(p, tp) {
  save_fig(p, paste0("fig4_gsea_bubble_", tp), width = 8, height = 6, subdir = "fig4")
})

save_fig(fig4_scatter, "fig4_innate_scatter_TX48", width = 7, height = 5, subdir = "fig4")
# ----

#+ FIGURE 7 - ALS organoid ageing timecourse
# Objects in data/fig7_als/ come from
# analysis/Anderson-Suthar_collab/data/r_objects/pre_checkpoint_objs
fig7_dir <- file.path(PROJ, "data", "fig7_als")

FIG7_NES_CUT  <- 1.5
FIG7_PVAL_CUT <- 0.05
FIG7_N_LABELS <- 25

# Fonts, pt
FIG7_AXIS_TITLE <- 20
FIG7_AXIS_TEXT  <- 16
FIG7_TITLE_SIZE <- 22
# Direction arrow labels, pt
FIG7_ARROW_TEXT <- 14
# Gene labels on the scatter, pt
FIG7_GENE_TEXT  <- 14

# Pathway labels truncated past this many characters
FIG7_LABEL_MAX <- 30

# Panel height per pathway row, relative to panel width
FIG7_ROW_ASPECT <- 0.085
# Bubble size in mm; -log10(padj) limits pooled across the four panels
FIG7_SIZE_RANGE  <- c(3, 9)
FIG7_SIZE_BREAKS_N <- 3

# Scatter x axis spans 1/FIG7_XMAX to FIG7_XMAX
FIG7_XMAX <- 4

# Scatter gene list
FIG7_INNATE_SETS <- c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_COMPLEMENT"
)

fig7_days <- c(d30 = "30", d50 = "50", d75 = "75", d120 = "120")

fig7_gsea <- setNames(lapply(names(fig7_days), function(d) {
  readRDS(file.path(fig7_dir, sprintf("fgsea_%s_ALS_vs_Ctrl(Batch).rds", d)))
}), names(fig7_days))

fig7_res <- readRDS(file.path(fig7_dir, "res_d120_ALS_vs_Ctrl(Batch).rds"))

# ---- Fig 7A-D: hallmark GSEA, C9orf72 vs Control
fig7_bubble_data <- function(gsea_res) {
  gsea_res |>
    dplyr::filter(abs(NES) > FIG7_NES_CUT, pval < FIG7_PVAL_CUT) |>
    dplyr::mutate(
      pathway_clean  = abbrev_pathway(pathway, FIG7_LABEL_MAX),
      direction      = ifelse(NES > 0, "Enriched", "Suppressed"),
      neg_log10_padj = -log10(padj)
    ) |>
    dplyr::arrange(NES)
}

fig7_bubble <- function(d, day_label, size_scale) {

  if (nrow(d) == 0) return(NULL)
  d$pathway_clean <- factor(d$pathway_clean, levels = d$pathway_clean)

  # Arrows sit in the space added below the first row; symmetric x limits put
  # one label in each half panel
  xm    <- max(abs(d$NES)) * 1.15
  y_arr <- -0.2
  y_lab <- -1.0

  ggplot(d, aes(NES, pathway_clean)) +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_clean),
                 colour = "grey60", linewidth = 0.4) +
    geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    annotate("segment", x = 0, xend = xm, y = y_arr, yend = y_arr,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             colour = "grey30", linewidth = 1.1) +
    annotate("segment", x = 0, xend = -xm, y = y_arr, yend = y_arr,
             arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
             colour = "grey30", linewidth = 1.1) +
    annotate("text", x = xm, y = y_lab, label = "Higher in C9orf72",
             hjust = 1, colour = "black", size = FIG7_ARROW_TEXT / .pt) +
    annotate("text", x = -xm, y = y_lab, label = "Higher in Control",
             hjust = 0, colour = "black", size = FIG7_ARROW_TEXT / .pt) +
    scale_fill_manual(values = c(Enriched = "#D73027", Suppressed = "#4575B4"),
                      guide = "none") +
    scale_size_continuous(range  = FIG7_SIZE_RANGE,
                          limits = size_scale$limits,
                          breaks = size_scale$breaks,
                          name   = "-log10(p-value)") +
    scale_x_continuous(limits = c(-xm, xm), expand = expansion(mult = 0.03)) +
    scale_y_discrete(expand = expansion(add = c(3.0, 0.6))) +
    labs(x = "Normalized Enrichment Score (NES)", y = NULL,
         title = sprintf("C9orf72 vs. Control (Day %s)", day_label)) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(),
          aspect.ratio = FIG7_ROW_ASPECT * (nrow(d) + 3),
          text         = element_text(colour = "black"),
          plot.title   = element_text(face = "bold", size = FIG7_TITLE_SIZE, colour = "black"),
          axis.title   = element_text(size = FIG7_AXIS_TITLE, colour = "black"),
          axis.text    = element_text(size = FIG7_AXIS_TEXT, colour = "black"),
          legend.text  = element_text(size = FIG7_AXIS_TEXT, colour = "black"),
          legend.title = element_text(size = FIG7_AXIS_TEXT, colour = "black"))
}

fig7_bubble_dfs <- lapply(fig7_gsea, fig7_bubble_data)

fig7_size <- bubble_size_scale(
  unlist(lapply(fig7_bubble_dfs, `[[`, "neg_log10_padj")),
  n = FIG7_SIZE_BREAKS_N
)

fig7_bubbles <- purrr::imap(fig7_bubble_dfs, function(d, day) {
  fig7_bubble(d, fig7_days[[day]], fig7_size)
})

fig7_bubbles$d30
fig7_bubbles$d50
fig7_bubbles$d75
fig7_bubbles$d120

# ---- Fig 7E: innate immunity genes at day 120
# Leading-edge genes of the innate sets, day 120 GSEA
fig7_targets <- fig7_gsea$d120 |>
  dplyr::filter(pathway %in% FIG7_INNATE_SETS) |>
  dplyr::pull(leadingEdge) |>
  unlist() |>
  unique()

fig7_volc <- fig7_res |>
  as.data.frame() |>
  rownames_to_column("gene") |>
  dplyr::mutate(
    FC        = 2^logFC,
    log10_p   = log10(P.Value),
    is_target = gene %in% fig7_targets,
    direction = dplyr::case_when(
      is_target & logFC > 0 ~ "Up",
      is_target & logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  )

# y limit set by the least significant highlighted gene
fig7_yfloor <- min(fig7_volc$log10_p[fig7_volc$is_target], na.rm = TRUE) * 1.2

# Band above the axis holding the direction arrows
fig7_ytop <- abs(fig7_yfloor) * 0.18

fig7_xarr <- FIG7_XMAX * 0.9

fig7_scatter <- ggplot(fig7_volc, aes(FC, log10_p)) +
  geom_point(data = function(d) dplyr::filter(d, !is_target),
             colour = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = function(d) dplyr::filter(d, is_target),
             aes(colour = direction), size = 1.8, alpha = 0.9) +
  # geom text size is mm, /.pt converts from pt
  geom_text_repel(data = function(d) dplyr::slice_min(dplyr::filter(d, is_target),
                                                      P.Value, n = FIG7_N_LABELS),
                  aes(label = gene), size = FIG7_GENE_TEXT / .pt, fontface = "italic",
                  max.overlaps = 20, segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.4,
                  ylim = c(fig7_yfloor, 0)) +
  geom_vline(xintercept = 1, linewidth = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  annotate("segment", x = 1, xend = fig7_xarr, y = fig7_ytop * 0.3, yend = fig7_ytop * 0.3,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           colour = "grey30", linewidth = 1.1) +
  annotate("segment", x = 1, xend = 1/fig7_xarr, y = fig7_ytop * 0.3, yend = fig7_ytop * 0.3,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           colour = "grey30", linewidth = 1.1) +
  annotate("text", x = fig7_xarr, y = fig7_ytop * 0.78, label = "Higher in C9orf72",
           hjust = 1, colour = "black", size = FIG7_ARROW_TEXT / .pt) +
  annotate("text", x = 1/fig7_xarr, y = fig7_ytop * 0.78, label = "Higher in Control",
           hjust = 0, colour = "black", size = FIG7_ARROW_TEXT / .pt) +
  scale_colour_manual(values = c(Up = "red", Down = "black"), guide = "none") +
  scale_x_continuous(trans = "log2", breaks = 2^(-5:5),
                     labels = c("-32", "-16", "-8", "-4", "-2", "1",
                                "+2", "+4", "+8", "+16", "+32")) +
  coord_cartesian(xlim = c(1/FIG7_XMAX, FIG7_XMAX), ylim = c(fig7_yfloor, fig7_ytop)) +
  labs(x = "FC: C9orf72 vs. Control", y = "log10(p-value)",
       title = "Innate Immunity Genes (Day 120)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        text       = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = FIG7_TITLE_SIZE, colour = "black"),
        axis.title = element_text(size = FIG7_AXIS_TITLE, colour = "black"),
        axis.text  = element_text(size = FIG7_AXIS_TEXT, colour = "black"))

fig7_scatter

# ---- Save Figure 7
purrr::iwalk(fig7_bubbles, function(p, d) {
  save_fig(p, paste0("fig7_gsea_bubble_", d), width = 8, height = 6, subdir = "fig7")
})

save_fig(fig7_scatter, "fig7_innate_scatter_d120", width = 7, height = 5, subdir = "fig7")
# ----

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
    # Legs from each sample to its group centroid
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

# ---- Interaction scatter ----
# Each interaction contrast as its two component responses: reference line's
# stimulus response on x, C9's on y. Off-diagonal points are the C9-specific
# response the contrast tests.
inter_reg <- contrast_registry |>
  dplyr::filter(type == "interaction", experiment == "spinal")

dir_cols <- c("Higher in C9" = "#D62728", "Lower in C9" = "#1F77B4", "NS" = "grey75")

inter_data <- purrr::pmap_dfr(
  list(inter_reg$name, as.character(inter_reg$stim), inter_reg$ref_line,
       inter_reg$label, inter_reg$min_n),
  function(nm, st, ref, lab, mn) {

    int <- tt$spinal[[nm]]
    c9  <- tt$spinal[[paste0("C9_", st, "_vs_Mock")]]
    ctl <- tt$spinal[[paste0(ref, "_", st, "_vs_Mock")]]

    tibble::tibble(
      contrast   = nm,
      label      = lab,
      stim       = st,
      ref_line   = ref,
      min_n      = mn,
      ensembl_id = int$ensembl_id,
      SYMBOL     = int$SYMBOL,
      x          = ctl$logFC[match(int$ensembl_id, ctl$ensembl_id)],
      y          = c9$logFC[match(int$ensembl_id, c9$ensembl_id)],
      int_lfc    = int$logFC,
      int_padj   = int$adj.P.Val
    ) |>
      dplyr::filter(!is.na(x), !is.na(y)) |>
      dplyr::mutate(
        sig = !is.na(int_padj) & int_padj < P_CUT & abs(int_lfc) > LFC_CUT,
        direction = dplyr::case_when(
          sig & int_lfc > 0 ~ "Higher in C9",
          sig & int_lfc < 0 ~ "Lower in C9",
          TRUE              ~ "NS"
        )
      )
  }
)

outdir <- ensure_dir(dir_fig, "interaction")

# One labelled panel per interaction contrast
for (nm in inter_reg$name) {

  df <- dplyr::filter(inter_data, contrast == nm)
  if (nrow(df) == 0) next

  lim <- range(c(df$x, df$y), finite = TRUE)
  r   <- stats::cor(df$x, df$y, use = "complete.obs")

  lab_df <- df |>
    dplyr::filter(sig, !is.na(SYMBOL), SYMBOL != "") |>
    dplyr::arrange(int_padj) |>
    head(15)

  caption <- if (df$min_n[1] == 1) {
    "Unreplicated (n = 1) - exploratory"
  } else {
    NULL
  }

  p <- ggplot(df, aes(x, y)) +
    geom_hline(yintercept = 0, colour = "grey88") +
    geom_vline(xintercept = 0, colour = "grey88") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
    geom_point(aes(colour = direction), size = 1.2, alpha = 0.6) +
    ggrepel::geom_text_repel(
      data = lab_df, aes(label = SYMBOL),
      size = 2.8, max.overlaps = 20, min.segment.length = 0,
      segment.colour = "grey60", segment.size = 0.3
    ) +
    scale_colour_manual(values = dir_cols, name = NULL) +
    coord_fixed(xlim = lim, ylim = lim) +
    labs(
      title    = df$label[1],
      subtitle = sprintf("r = %.2f   |   %d genes differ (adj.P.Val < %s, |logFC| > %s)",
                         r, sum(df$sig), P_CUT, LFC_CUT),
      caption  = caption,
      x = paste0(df$ref_line[1], ": ", df$stim[1], " vs Mock (log2FC)"),
      y = paste0("C9: ", df$stim[1], " vs Mock (log2FC)")
    ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.caption = element_text(colour = "#B22222", hjust = 0))

  ggsave(file.path(outdir, paste0("interaction_", nm, ".png")),
         p, width = 6.5, height = 6.5, dpi = 300)
}

# Overview grid: reference line x stimulus, no gene labels
inter_overview <- inter_data |>
  dplyr::mutate(
    stim     = factor(stim, levels = STIM_LEVELS),
    ref_line = factor(ref_line, levels = LINE_LEVELS)
  )

lim <- range(c(inter_overview$x, inter_overview$y), finite = TRUE)

p_all <- ggplot(inter_overview, aes(x, y)) +
  geom_hline(yintercept = 0, colour = "grey88") +
  geom_vline(xintercept = 0, colour = "grey88") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey35") +
  geom_point(aes(colour = direction), size = 0.9, alpha = 0.55) +
  facet_grid(ref_line ~ stim) +
  scale_colour_manual(values = dir_cols, name = NULL) +
  coord_fixed(xlim = lim, ylim = lim) +
  labs(
    title = "C9 stimulus response vs control line response",
    subtitle = "Off-diagonal genes drive the interaction contrasts",
    x = "Control line: stimulus vs Mock (log2FC)",
    y = "C9: stimulus vs Mock (log2FC)"
  ) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = NA))

ggsave(file.path(outdir, "interaction_overview.png"), p_all,
       width = 11, height = 8, dpi = 300)

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
# EDIT: pick a hallmark set from script 4's GSEA results
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
