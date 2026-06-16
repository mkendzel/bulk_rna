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
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD.rds")

# ---- Rename ----
res_shrunk_list <- purrr::imap(res_shrunk_list, ~{
  
  df <- as.data.frame(.x)
  df$ensembl_id <- rownames(df)
  
  n_before <- nrow(df)
  
  map_df <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys     = df$ensembl_id,
    keytype  = "ENSEMBL",
    columns  = "SYMBOL"
  )
  
  map_df <- as.data.frame(map_df)
  map_df <- map_df[!duplicated(map_df$ENSEMBL), ]
  
  df <- merge(df, map_df,
              by.x = "ensembl_id",
              by.y = "ENSEMBL",
              all.x = TRUE,
              sort = FALSE)
  
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "", ]
  
  n_after <- nrow(df)
  removed <- n_before - n_after
  
  message(.y, ": removed ", removed, " genes (", n_before, " → ", n_after, ")")
  
  rownames(df) <- df$ensembl_id
  df
})

# ---- order ----
ranked_list <- res_shrunk_list %>%
  purrr::map(~{
    
    df <- as.data.frame(.x)
    df <- df[!is.na(df$SYMBOL), ]
    
    ranks_logfc <- df$log2FoldChange
    names(ranks_logfc) <- df$SYMBOL
    
    ranks_logfc <- ranks_logfc[is.finite(ranks_logfc)]
    ranks_logfc <- ranks_logfc[!duplicated(names(ranks_logfc))]
    ranks_logfc <- sort(ranks_logfc, decreasing = TRUE)
    
    ranks_logfc
  })


gsea_results <- map(
  ranked_list,
  ~ fgsea(
    pathways = hallmark_sets,
    stats    = .x
  )
)

# ---- GSEA Results ----

gsea_tbl_list <- imap(gsea_results, ~{
  as_tibble(.x) %>%
    arrange(padj) %>%
    mutate(contrast = .y, .before = 1)
})

# ---- Bubble graph ----
# Parse names like: condition_ROV_12_vs_Mock_0
parse_contrast <- function(comparison) {
  x <- sub("^condition_", "", comparison)
  parts <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
  treat <- parts[1]
  ctrl  <- if (length(parts) > 1) parts[2] else NA_character_
  
  tr <- strsplit(treat, "_", fixed = TRUE)[[1]]
  treatment <- tr[1]
  time      <- suppressWarnings(as.integer(tr[2]))
  
  list(
    treat = treat,
    ctrl = ctrl,
    treatment = treatment,
    time = time
  )
}

# Combine all contrasts into one table
df_plot <- purrr::imap_dfr(gsea_results, ~{
  as_tibble(.x) %>%
    mutate(comparison = .y, .before = 1)
})

# Parse labels for x axis
ps <- lapply(df_plot$comparison, parse_contrast)
df_plot$treat     <- vapply(ps, `[[`, character(1), "treat")
df_plot$ctrl       <- vapply(ps, `[[`, character(1), "ctrl")
df_plot$treatment  <- vapply(ps, `[[`, character(1), "treatment")
df_plot$time       <- vapply(ps, `[[`, integer(1),   "time")

# Clean pathway labels (Hallmark)
df_plot$pathway <- sub("^HALLMARK_", "", df_plot$pathway)

# Filter to significant pathways only
df_plot <- df_plot %>%
  filter(!is.na(padj), padj <= 0.05, !is.na(NES), !is.na(pathway))

# Order x-axis by time, then treatment (within time)
x_levels <- df_plot %>%
  distinct(treat, treatment, time) %>%
  arrange(time, treatment) %>%
  pull(treat)

df_plot$treat <- factor(df_plot$treat, levels = x_levels)

# Order y-axis by mean NES
df_plot$pathway <- with(df_plot, reorder(pathway, NES, mean))

# Bubble plot
GSEA_bubble <- ggplot(
  df_plot,
  aes(x = treat, y = pathway, size = -log10(padj), fill = NES)
) +
  geom_point(shape = 21, color = "black") +
  scale_size(name = "-log10(padj)", range = c(4, 10)) +
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

GSEA_bubble

# ---- Other organoids ----
