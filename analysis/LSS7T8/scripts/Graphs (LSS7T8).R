# ---- Library and helper functions ----
library(tidyverse)
library(fgsea)
library(ggrepel)
library(patchwork)

# source helper functions from R/
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

ckpt_dir <- "projects/LSS7T8/data/r_objects"
fig_dir  <- "projects/LSS7T8/figures"
pathways <- gmtPathways("genesets/mh.all.v2026.1.Mm.symbols.gmt")

# Significance cutoffs shared by the volcano guides, point colours, and subtitles
padj_cutoff <- 0.05
lfc_cutoff  <- 1

# Save a ggplot to fig_dir as PDF, dimensions in inches.
# cairo_pdf because labels use non-ASCII glyphs (alpha, gamma, arrows).
save_fig <- function(plot, filename, width, height) {
  ggsave(file.path(fig_dir, filename), plot,
         width = width, height = height, device = cairo_pdf)
}

# limma-voom DE results: named list, one topTable per contrast
res_liver  <- load_checkpoint("res_voom_liver_annotated", dir = ckpt_dir)
res_spleen <- load_checkpoint("res_voom_spleen_annotated", dir = ckpt_dir)

# Draft watermark layer, centred on the panel. Add last so it draws on top.
watermark <- function(label = "Not for publication") {
  annotation_custom(
    grid::textGrob(
      label,
      rot = 30,
      gp = grid::gpar(fontsize = 24, col = "grey50", alpha = 0.18, fontface = "bold")
    )
  )
}

# X-axis label with arrows marking which group each side is higher in.
direction_axis_label <- function(value_label, group_neg, group_pos) {
  paste0(value_label, "\n← Higher in ", group_neg, "   |   Higher in ", group_pos, " →")
}

# Keep the top-ranked row per cell of an n_bins x n_bins grid over the plotted
# range. Crowded regions get thinned to one label, isolated points always survive.
thin_labels <- function(df, x, y, rank_by, n_bins = 8) {
  if (nrow(df) < 2) return(df)

  df %>%
    mutate(
      .x_bin = cut({{ x }}, breaks = n_bins, labels = FALSE),
      .y_bin = cut({{ y }}, breaks = n_bins, labels = FALSE)
    ) %>%
    group_by(.x_bin, .y_bin) %>%
    slice_max({{ rank_by }}, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(-.x_bin, -.y_bin)
}

# Candidate labels for a volcano: significant genes, one row per symbol, thinned
# on the grid. Rank favours genes that are both large-effect and confident.
label_set <- function(volcano_df, n_bins) {
  volcano_df %>%
    filter(sig) %>%
    mutate(label_rank = neg_log10_padj * abs(logFC)) %>%
    arrange(desc(label_rank)) %>%
    distinct(SYMBOL, .keep_all = TRUE) %>%
    thin_labels(logFC, neg_log10_padj, label_rank, n_bins = n_bins)
}

# Volcano of one topTable. Grey = not in pathway_name, blue = in it,
# red = in it and past the dashed cutoffs. Labels are thinned. y capped at 50.
make_ifn_volcano <- function(res_df, pathway_name, set_label, contrast_label, tissue,
                             group_neg, group_pos, n_bins = 8) {
  target_genes <- pathways[[pathway_name]]

  volcano_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(adj.P.Val), 50),
      is_target = SYMBOL %in% target_genes,
      sig = is_target & adj.P.Val < padj_cutoff & abs(logFC) > lfc_cutoff
    )

  n_detected <- n_distinct(volcano_df$SYMBOL[volcano_df$is_target])
  n_sig      <- n_distinct(volcano_df$SYMBOL[volcano_df$sig])
  x_limit    <- max(abs(volcano_df$logFC)) * 1.05

  ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(data = filter(volcano_df, !is_target),
               color = "grey75", size = 0.5, alpha = 0.4) +
    geom_point(data = filter(volcano_df, is_target & !sig),
               color = "steelblue", size = 1.5, alpha = 0.7) +
    geom_point(data = filter(volcano_df, sig),
               color = "firebrick", size = 2, alpha = 0.9) +
    geom_text_repel(
      data = label_set(volcano_df, n_bins),
      aes(label = SYMBOL),
      size = 3, max.overlaps = Inf, segment.color = "grey40",
      color = "firebrick", seed = 42
    ) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      x = direction_axis_label("log2 Fold Change", group_neg, group_pos),
      y = "-log10(adj.P.Val)",
      title = paste0(contrast_label, " (", tissue, ")"),
      subtitle = paste0(set_label, " — ", n_sig, "/", n_detected, " set genes significant")
    ) +
    watermark() +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

# Volcano of every gene in a topTable, coloured by direction. Same cutoffs,
# labels, and y cap as make_ifn_volcano.
make_volcano <- function(res_df, contrast_label, tissue, group_neg, group_pos, n_bins = 8) {
  volcano_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      neg_log10_padj = pmin(-log10(adj.P.Val), 50),
      sig = adj.P.Val < padj_cutoff & abs(logFC) > lfc_cutoff,
      direction = factor(
        case_when(sig & logFC > 0 ~ group_pos,
                  sig & logFC < 0 ~ group_neg,
                  TRUE ~ "Not significant"),
        levels = c(group_neg, group_pos, "Not significant")
      )
    )

  n_up    <- n_distinct(volcano_df$SYMBOL[volcano_df$sig & volcano_df$logFC > 0])
  n_down  <- n_distinct(volcano_df$SYMBOL[volcano_df$sig & volcano_df$logFC < 0])
  x_limit <- max(abs(volcano_df$logFC)) * 1.05

  color_values <- setNames(c("#4575B4", "#D73027", "grey75"),
                           c(group_neg, group_pos, "Not significant"))

  ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(data = filter(volcano_df, !sig),
               aes(color = direction), size = 0.5, alpha = 0.4) +
    geom_point(data = filter(volcano_df, sig),
               aes(color = direction), size = 1.5, alpha = 0.8) +
    geom_text_repel(
      data = label_set(volcano_df, n_bins),
      aes(label = SYMBOL, color = direction),
      size = 3, max.overlaps = Inf, segment.color = "grey40",
      show.legend = FALSE, seed = 42
    ) +
    scale_color_manual(values = color_values, drop = FALSE) +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    labs(
      x = direction_axis_label("log2 Fold Change", group_neg, group_pos),
      y = "-log10(adj.P.Val)",
      color = "Higher in",
      title = paste0(contrast_label, " (", tissue, ")"),
      subtitle = paste0("All genes — ", n_up, " higher in ", group_pos, ", ",
                        n_down, " higher in ", group_neg,
                        " (adj.P < ", padj_cutoff, ", |logFC| > ", lfc_cutoff, ")")
    ) +
    watermark() +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12),
          plot.subtitle = element_text(size = 9, color = "grey30"))
}

# Lollipop of fgsea pathways with |NES| > 1.5 and pval < 0.05, ordered by NES.
# Point size is -log10(padj). NES < 0 is higher in group_neg, NES > 0 in group_pos.
make_lollipop <- function(gsea_df, plot_title, group_neg, group_pos) {
  sig_paths <- gsea_df %>%
    filter(abs(NES) > 1.5, pval < 0.05) %>%
    mutate(
      pathway_clean = str_remove(pathway, "^HALLMARK_") %>%
        str_replace_all("_", " ") %>%
        str_to_title(),
      direction = factor(ifelse(NES > 0, group_pos, group_neg),
                         levels = c(group_neg, group_pos)),
      neg_log10_padj = -log10(padj)
    ) %>%
    arrange(NES)

  sig_paths$pathway_clean <- factor(sig_paths$pathway_clean, levels = sig_paths$pathway_clean)

  fill_values <- setNames(c("#4575B4", "#D73027"), c(group_neg, group_pos))

  x_lab <- direction_axis_label("Normalized Enrichment Score (NES)", group_neg, group_pos)

  ggplot(sig_paths, aes(x = NES, y = pathway_clean)) +
    geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean),
                 color = "grey60", linewidth = 0.4) +
    geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
    scale_fill_manual(values = fill_values, drop = FALSE) +
    scale_size_continuous(range = c(2, 7), name = "-log10(padj)") +
    scale_x_continuous(limits = max(abs(sig_paths$NES)) * c(-1.05, 1.05)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    labs(x = x_lab, y = NULL, fill = "Higher in", title = plot_title) +
    watermark() +
    theme_minimal(base_size = 11) +
    theme(panel.grid.major.y = element_blank(), legend.position = "right",
          plot.title = element_text(face = "bold", size = 12))
}


# ================================ LIVER ================================

# ---- CL vs Control ----

gsea_results <- load_checkpoint("fgsea_cl_vs_mock_liver", dir = ckpt_dir)
res_results  <- res_liver[["treatmentcl"]]

tissue    <- "Liver"
contrast  <- "CL vs. Control"
group_neg <- "Control"
group_pos <- "CL"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_all <- make_volcano(res_results, contrast, tissue, group_neg, group_pos)

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

p_volcano_alpha

save_fig(p_lollipop,      "lollipop_cl_vs_mock_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_all,   "volcano_all_cl_vs_mock_liver.pdf",  width = 7, height = 5)
save_fig(p_volcano_alpha, "volcano_ifna_cl_vs_mock_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_cl_vs_mock_liver.pdf", width = 7, height = 5)


# ---- DOPC vs Control ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_mock_liver", dir = ckpt_dir)
res_results  <- res_liver[["treatmentdopc"]]

tissue    <- "Liver"
contrast  <- "DOPC vs. Control"
group_neg <- "Control"
group_pos <- "DOPC"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_all <- make_volcano(res_results, contrast, tissue, group_neg, group_pos)

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

save_fig(p_lollipop,      "lollipop_dopc_vs_mock_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_all,   "volcano_all_dopc_vs_mock_liver.pdf",  width = 7, height = 5)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_mock_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_mock_liver.pdf", width = 7, height = 5)


# ---- DOPC vs CL ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_cl_liver", dir = ckpt_dir)
res_results  <- res_liver[["treatment_dopc_vs_cl"]]

tissue    <- "Liver"
contrast  <- "DOPC vs. CL"
group_neg <- "CL"
group_pos <- "DOPC"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

save_fig(p_lollipop,      "lollipop_dopc_vs_cl_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_cl_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_cl_liver.pdf", width = 7, height = 5)


# ================================ SPLEEN ================================

# ---- CL vs Control ----

gsea_results <- load_checkpoint("fgsea_cl_vs_mock_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatmentcl"]]

tissue    <- "Spleen"
contrast  <- "CL vs. Control"
group_neg <- "Control"
group_pos <- "CL"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_all <- make_volcano(res_results, contrast, tissue, group_neg, group_pos)

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

p_volcano_alpha

save_fig(p_lollipop,      "lollipop_cl_vs_mock_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_all,   "volcano_all_cl_vs_mock_spleen.pdf",  width = 7, height = 5)
save_fig(p_volcano_alpha, "volcano_ifna_cl_vs_mock_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_cl_vs_mock_spleen.pdf", width = 7, height = 5)


# ---- DOPC vs Control ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_mock_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatmentdopc"]]

tissue    <- "Spleen"
contrast  <- "DOPC vs. Control"
group_neg <- "Control"
group_pos <- "DOPC"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_all <- make_volcano(res_results, contrast, tissue, group_neg, group_pos)

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

save_fig(p_lollipop,      "lollipop_dopc_vs_mock_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_all,   "volcano_all_dopc_vs_mock_spleen.pdf",  width = 7, height = 5)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_mock_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_mock_spleen.pdf", width = 7, height = 5)


# ---- DOPC vs CL ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_cl_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatment_dopc_vs_cl"]]

tissue    <- "Spleen"
contrast  <- "DOPC vs. CL"
group_neg <- "CL"
group_pos <- "DOPC"

p_lollipop <- make_lollipop(gsea_results, paste0(contrast, " (", tissue, ")"),
                            group_neg, group_pos)

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response", contrast, tissue, group_neg, group_pos
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response", contrast, tissue, group_neg, group_pos
)

save_fig(p_lollipop,      "lollipop_dopc_vs_cl_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_cl_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_cl_spleen.pdf", width = 7, height = 5)
