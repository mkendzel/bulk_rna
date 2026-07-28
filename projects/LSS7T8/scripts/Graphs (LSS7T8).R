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

# Volcano of one topTable. Grey = not in pathway_name, blue = in it,
# red = in it and past the dashed cutoffs (adj.P.Val < 0.05, |logFC| > 1).
# Red genes are labelled. y capped at 50.
make_ifn_volcano <- function(res_df, pathway_name, plot_title) {
  target_genes <- pathways[[pathway_name]]

  volcano_df <- res_df %>%
    filter(!is.na(SYMBOL), !is.na(logFC), !is.na(adj.P.Val)) %>%
    mutate(
      neg_log10_padj = -log10(adj.P.Val),
      sig = adj.P.Val < 0.05 & abs(logFC) > 1,
      is_target = SYMBOL %in% target_genes,
      neg_log10_padj = pmin(neg_log10_padj, 50)
    )

  ggplot(volcano_df, aes(x = logFC, y = neg_log10_padj)) +
    geom_point(data = filter(volcano_df, !is_target),
               color = "grey75", size = 0.5, alpha = 0.4) +
    geom_point(data = filter(volcano_df, is_target & !sig),
               color = "steelblue", size = 1.5, alpha = 0.7) +
    geom_point(data = filter(volcano_df, is_target & sig),
               color = "firebrick", size = 2, alpha = 0.9) +
    geom_text_repel(
      data = filter(volcano_df, is_target & sig),
      aes(label = SYMBOL),
      size = 3, max.overlaps = 25, segment.color = "grey40",
      color = "firebrick"
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    labs(x = "log2 Fold Change", y = "-log10(adj.P.Val)", title = plot_title) +
    watermark() +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold", size = 12))
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

  x_lab <- paste0(
    "Normalized Enrichment Score (NES)\n",
    "← Higher in ", group_neg, "   |   Higher in ", group_pos, " →"
  )

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

p_lollipop <- make_lollipop(gsea_results, "CL vs. Control (Liver)", "Control", "CL")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

p_volcano_alpha

save_fig(p_lollipop,      "lollipop_cl_vs_mock_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_cl_vs_mock_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_cl_vs_mock_liver.pdf", width = 7, height = 5)


# ---- DOPC vs Control ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_mock_liver", dir = ckpt_dir)
res_results  <- res_liver[["treatmentdopc"]]

p_lollipop <- make_lollipop(gsea_results, "DOPC vs. Control (Liver)", "Control", "DOPC")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

save_fig(p_lollipop,      "lollipop_dopc_vs_mock_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_mock_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_mock_liver.pdf", width = 7, height = 5)


# ---- DOPC vs CL ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_cl_liver", dir = ckpt_dir)
res_results  <- res_liver[["treatment_dopc_vs_cl"]]

p_lollipop <- make_lollipop(gsea_results, "DOPC vs. CL (Liver)", "CL", "DOPC")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

save_fig(p_lollipop,      "lollipop_dopc_vs_cl_liver.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_cl_liver.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_cl_liver.pdf", width = 7, height = 5)


# ================================ SPLEEN ================================

# ---- CL vs Control ----

gsea_results <- load_checkpoint("fgsea_cl_vs_mock_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatmentcl"]]

p_lollipop <- make_lollipop(gsea_results, "CL vs. Control (Spleen)", "Control", "CL")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

p_volcano_alpha

save_fig(p_lollipop,      "lollipop_cl_vs_mock_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_cl_vs_mock_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_cl_vs_mock_spleen.pdf", width = 7, height = 5)


# ---- DOPC vs Control ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_mock_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatmentdopc"]]

p_lollipop <- make_lollipop(gsea_results, "DOPC vs. Control (Spleen)", "Control", "DOPC")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

save_fig(p_lollipop,      "lollipop_dopc_vs_mock_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_mock_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_mock_spleen.pdf", width = 7, height = 5)


# ---- DOPC vs CL ----

gsea_results <- load_checkpoint("fgsea_dopc_vs_cl_spleen", dir = ckpt_dir)
res_results  <- res_spleen[["treatment_dopc_vs_cl"]]

p_lollipop <- make_lollipop(gsea_results, "DOPC vs. CL (Spleen)", "CL", "DOPC")

p_lollipop

p_volcano_alpha <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "IFN-α Response Genes"
)

p_volcano_gamma <- make_ifn_volcano(
  res_results,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "IFN-γ Response Genes"
)

save_fig(p_lollipop,      "lollipop_dopc_vs_cl_spleen.pdf",     width = 8, height = 6)
save_fig(p_volcano_alpha, "volcano_ifna_dopc_vs_cl_spleen.pdf", width = 7, height = 5)
save_fig(p_volcano_gamma, "volcano_ifng_dopc_vs_cl_spleen.pdf", width = 7, height = 5)
