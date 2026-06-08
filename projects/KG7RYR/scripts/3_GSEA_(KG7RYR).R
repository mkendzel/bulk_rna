# ---- Libraries ----
library(Matrix)
library(DESeq2)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)
library(pheatmap)
library(grid)
library(knitr)
library(openxlsx)
library(purrr)
library(tibble)
library(ggplot2)

# ---- Dataset ----
gmt_path <- "data/geneset/mh.all.v2026.1.Mm.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

# liver Results
res_shrunk_liver <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_liver_qcmin10_annotated.rds")
res_wald_liver <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_liver_qcmin10_annotated.rds")

# Spleen Results
res_wald_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")
res_shrunk_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")

# ---- GSEA for liver ----
ranked_list_liver <- res_shrunk_liver %>%
  purrr::map(~{
    
    df <- as.data.frame(.x)
    df <- df[!is.na(df$symbol), ]
    
    ranks_logfc <- df$log2FoldChange
    names(ranks_logfc) <- df$symbol
    
    ranks_logfc <- ranks_logfc[is.finite(ranks_logfc)]
    ranks_logfc <- ranks_logfc[!duplicated(names(ranks_logfc))]
    ranks_logfc <- sort(ranks_logfc, decreasing = TRUE)
    
    ranks_logfc
  })

ranked_list_liver <- res_shrunk_liver %>%
  purrr::map(~{
    
    df <- as.data.frame(.x)
    df <- df[!is.na(df$symbol), ]
    
    ranks <- sign(df$log2FoldChange) * -log10(df$pvalue)
    names(ranks) <- df$symbol
    
    ranks <- ranks[is.finite(ranks)]
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)
    
    ranks
  })

gsea_results_liver <- map(
  ranked_list_liver,
  ~ fgsea(
    pathways = hallmark_sets,
    stats    = .x
  )
)

saveRDS(gsea_results_liver, "data/KG7RYR/r_objects/gsea_results_liver(signLog_padj).rds")
# Parse names like
parse_contrast_liver <- function(comparison) {
  x <- sub("^treatment_", "", comparison)
  parts <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
  treat <- parts[1]
  ctrl  <- if (length(parts) > 1) parts[2] else NA_character_
  
  list(
    treat = treat,
    ctrl = ctrl
  )
}

df_plot_liver <- purrr::imap_dfr(gsea_results_liver, ~{
  as_tibble(.x) %>%
    mutate(comparison = .y, .before = 1)
})

# Parse labels for x axis
ps <- lapply(df_plot_liver$comparison, parse_contrast_liver)
df_plot_liver$treat <- vapply(ps, `[[`, character(1), "treat")
df_plot_liver$ctrl  <- vapply(ps, `[[`, character(1), "ctrl")

# Clean pathway labels (Hallmark)
df_plot_liver$pathway <- sub("^HALLMARK_", "", df_plot_liver$pathway)

# Filter to significant pathways only
df_plot_liver <- df_plot_liver %>%
  filter(!is.na(padj), padj <= 0.05, !is.na(NES), !is.na(pathway))

# Order x-axis
x_levels_liver <- df_plot_liver %>%
  distinct(treat) %>%
  arrange(treat) %>%
  pull(treat)

df_plot_liver$treat <- factor(df_plot_liver$treat, levels = x_levels_liver)

#Order y-axis
df_plot_liver$pathway <- with(df_plot_liver, reorder(pathway, NES, mean))

df_plot_liver$neg_log10_padj <- pmin(-log10(df_plot_liver$padj), 10)

GSEA_bubble_liver <- ggplot(
  df_plot_liver,
  aes(x = treat, y = pathway, size = neg_log10_padj, fill = NES)
) +
  geom_point(shape = 21, color = "black") +
  scale_size(name = "-log10(padj)", range = c(4, 10)) +
  ggtitle("Liver") +
  ylab("Gene Set") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 0),
    axis.title.x = element_blank()
  ) +
  scale_fill_distiller(palette = "Spectral")
GSEA_bubble_liver

# ---- Spleen GSEA ----
ranked_list_spleen <- res_shrunk_spleen %>%
  purrr::map(~{
    
    df <- as.data.frame(.x)
    df <- df[!is.na(df$symbol), ]
    
    ranks <- sign(df$log2FoldChange) * -log10(df$pvalue)
    names(ranks) <- df$symbol
    
    ranks <- ranks[is.finite(ranks)]
    ranks <- ranks[!duplicated(names(ranks))]
    ranks <- sort(ranks, decreasing = TRUE)
    
    ranks
  })

gsea_results_spleen <- map(
  ranked_list_spleen,
  ~ fgsea(
    pathways = hallmark_sets,
    stats    = .x
  )
)
saveRDS(gsea_results_spleen, "data/KG7RYR/r_objects/gsea_results_spleen(signLog_padj).rds")

parse_contrast_spleen <- function(comparison) {
  x <- sub("^treatment_", "", comparison)
  parts <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
  treat <- parts[1]
  ctrl  <- if (length(parts) > 1) parts[2] else NA_character_
  
  list(
    treat = treat,
    ctrl = ctrl
  )
}

df_plot_spleen <- purrr::imap_dfr(gsea_results_spleen, ~{
  as_tibble(.x) %>%
    mutate(comparison = .y, .before = 1)
})

ps <- lapply(df_plot_spleen$comparison, parse_contrast_spleen)
df_plot_spleen$treat <- vapply(ps, `[[`, character(1), "treat")
df_plot_spleen$ctrl  <- vapply(ps, `[[`, character(1), "ctrl")

df_plot_spleen$pathway <- sub("^HALLMARK_", "", df_plot_spleen$pathway)

df_plot_spleen <- df_plot_spleen %>%
  filter(!is.na(padj), padj <= 0.05, !is.na(NES), !is.na(pathway))

x_levels_spleen <- df_plot_spleen %>%
  distinct(treat) %>%
  arrange(treat) %>%
  pull(treat)
df_plot_spleen$treat <- factor(df_plot_spleen$treat, levels = x_levels_spleen)

df_plot_spleen$pathway <- with(df_plot_spleen, reorder(pathway, NES, mean))
df_plot_spleen$neg_log10_padj <- pmin(-log10(df_plot_spleen$padj), 10)

GSEA_bubble_spleen <- ggplot(
  df_plot_spleen,
  aes(x = treat, y = pathway, size = neg_log10_padj, fill = NES)
) +
  geom_point(shape = 21, color = "black") +
  scale_size(name = "-log10(padj)", range = c(4, 10)) +
  ggtitle("Spleen") +
  ylab("Gene Set") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 0),
    axis.title.x = element_blank()
  ) +
  scale_fill_distiller(palette = "Spectral")
GSEA_bubble_spleen
