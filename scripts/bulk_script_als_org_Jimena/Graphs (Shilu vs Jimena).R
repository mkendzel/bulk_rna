# ---- Data Load ----
# ---- GSEA graphs ----
library(ggplot2)
library(dplyr)

# Helper to prep fgsea results for plotting
prep_gsea_plot <- function(gsea_res, padj_cutoff = 0.05) {
  gsea_res %>%
    as.data.frame() %>%
    filter(padj < padj_cutoff) %>%
    mutate(
      Direction = ifelse(NES > 0, "Up", "Down"),
      neg_log10_padj = -log10(padj),
      # Clean up pathway names for display
      pathway = gsub("^HALLMARK_", "", pathway),
      pathway = gsub("_", " ", pathway)
    ) %>%
    arrange(NES)
}

# Plotting function
plot_gsea <- function(gsea_res, title, padj_cutoff = 0.05) {
  df <- prep_gsea_plot(gsea_res, padj_cutoff)
  
  if (nrow(df) == 0) {
    message("No significant pathways at padj < ", padj_cutoff)
    return(NULL)
  }
  
  df$pathway <- factor(df$pathway, levels = df$pathway)
  
  ggplot(df, aes(x = NES, y = pathway, color = Direction)) +
    geom_segment(aes(x = 0, xend = NES, y = pathway, yend = pathway),
                 color = "grey60", linewidth = 0.7) +
    geom_point(aes(size = neg_log10_padj)) +
    scale_color_manual(values = c("Down" = "steelblue", "Up" = "firebrick")) +
    labs(
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      size = "-log10(adj p)",
      title = title
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(face = "bold", size = 12)
    )
}

# Generate plots
p_d30  <- plot_gsea(gsea_d30,  "GSEA - Up & Down Pathways (d30)")
p_d50  <- plot_gsea(gsea_d50,  "GSEA - Up & Down Pathways (d50)")
p_d75  <- plot_gsea(gsea_d75,  "GSEA - Up & Down Pathways (d75)")
p_d120 <- plot_gsea(gsea_d120, "GSEA - Up & Down Pathways (d120)")

# Save
ggsave("figures/Shilu's vs Jimena/GSEA_d30.pdf", p_d30, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d50.pdf", p_d50, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d75.pdf", p_d75, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d120.pdf", p_d120, width = 8, height = 6)

# ---- Network activity of GSEA ----

library(ggplot2)
library(ggrepel)
library(dplyr)
library(enrichplot)
library(clusterProfiler)

# Data

cp_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d30_ALS_vs_Ctrl.rds")
cp_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d50_ALS_vs_Ctrl.rds")
cp_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d75_ALS_vs_Ctrl.rds")
cp_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d120_ALS_vs_Ctrl.rds")

res_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d30_ALS_vs_Ctrl.rds")
res_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d50_ALS_vs_Ctrl.rds")
res_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d75_ALS_vs_Ctrl.rds")
res_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d120_ALS_vs_Ctrl.rds")

#plot function

plot_cnet <- function(cp_res, tt, title, show = 5,
                      filename = NULL, width = 16, height = 14) {
  
  p <- cnetplot(cp_res, showCategory = show)
  
  gene_data <- p$data[!p$data$.isCategory, ]
  path_data <- p$data[p$data$.isCategory, ]
  
  lfc <- setNames(tt$logFC, rownames(tt))
  lfc <- lfc[!duplicated(names(lfc))]
  gene_data$lfc <- lfc[gene_data$name]
  
  top_paths <- cp_res@result %>% arrange(p.adjust) %>% slice_head(n = show)
  pathway_nes <- setNames(top_paths$NES, top_paths$ID)
  path_data$direction <- ifelse(pathway_nes[path_data$name] > 0, "Up", "Down")
  
  path_data$label <- gsub("^HALLMARK_", "", path_data$name)
  path_data$label <- gsub("_", " ", path_data$label)
  
  all_labels <- rbind(
    data.frame(x = gene_data$x, y = gene_data$y, label = gene_data$name,
               is_pathway = FALSE, stringsAsFactors = FALSE),
    data.frame(x = path_data$x, y = path_data$y, label = path_data$label,
               is_pathway = TRUE, stringsAsFactors = FALSE)
  )
  
  p$layers[[4]] <- NULL
  p$layers[[3]] <- NULL
  p$layers[[2]] <- NULL
  
  p <- p +
    geom_point(
      data = gene_data,
      aes(x = x, y = y, color = lfc),
      size = 3,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = path_data,
      aes(x = x, y = y),
      size = 8,
      fill = ifelse(path_data$direction == "Up", "#D55E00", "#0072B2"),
      color = "black",
      shape = 21,
      stroke = 1.2,
      inherit.aes = FALSE
    ) +
    geom_text_repel(
      data = all_labels,
      aes(x = x, y = y, label = label),
      fontface = ifelse(all_labels$is_pathway, "bold", "plain"),
      size = ifelse(all_labels$is_pathway, 4, 3),
      max.overlaps = Inf,
      box.padding = 0.4,
      point.padding = 0.3,
      min.segment.length = 0.2,
      segment.color = "grey50",
      inherit.aes = FALSE
    ) +
    scale_color_gradient2(low = "#0072B2", mid = "white", high = "#D55E00",
                          midpoint = 0, name = "logFC") +
    ggtitle(title)
  
  if (!is.null(filename)) {
    ggsave(filename, p, width = width, height = height)
  }
  
  return(p)
}

# Run and save all four
plot_cnet(cp_d30,  res_d30,  "GSEA Network (d30 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d30.pdf")
plot_cnet(cp_d50,  res_d50,  "GSEA Network (d50 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d50.pdf")
plot_cnet(cp_d75,  res_d75,  "GSEA Network (d75 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d75.pdf")
plot_cnet(cp_d120, res_d120, "GSEA Network (d120 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d120.pdf")