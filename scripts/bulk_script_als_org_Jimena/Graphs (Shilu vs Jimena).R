# ---- GSEA graphs ----
library(ggplot2)
library(dplyr)

gsea_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d30_ALS_vs_Ctrl(Batch).rds")
gsea_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d50_ALS_vs_Ctrl(Batch).rds")
gsea_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d75_ALS_vs_Ctrl(Batch).rds")
gsea_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d120_ALS_vs_Ctrl(Batch).rds")
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
ggsave("figures/Shilu's vs Jimena/GSEA_d30(Batch).pdf", p_d30, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d50(Batch).pdf", p_d50, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d75(Batch).pdf", p_d75, width = 8, height = 6)
ggsave("figures/Shilu's vs Jimena/GSEA_d120(Batch).pdf", p_d120, width = 8, height = 6)

# ---- Network activity of GSEA ----

library(ggplot2)
library(ggrepel)
library(dplyr)
library(enrichplot)
library(clusterProfiler)

# Data

cp_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d30_ALS_vs_Ctrl(Batch).rds")
cp_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d50_ALS_vs_Ctrl(Batch).rds")
cp_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d75_ALS_vs_Ctrl(Batch).rds")
cp_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/cp_gsea_d120_ALS_vs_Ctrl(Batch).rds")

res_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d30_ALS_vs_Ctrl(Batch).rds")
res_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d50_ALS_vs_Ctrl(Batch).rds")
res_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d75_ALS_vs_Ctrl(Batch).rds")
res_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d120_ALS_vs_Ctrl(Batch).rds")

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
          filename = "figures/Shilu's vs Jimena/cnet_d30(Batch).pdf")
plot_cnet(cp_d50,  res_d50,  "GSEA Network (d50 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d50(Batch).pdf")
plot_cnet(cp_d75,  res_d75,  "GSEA Network (d75 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d75(Batch).pdf")
plot_cnet(cp_d120, res_d120, "GSEA Network (d120 - ALS vs Control)",
          filename = "figures/Shilu's vs Jimena/cnet_d120(Batch).pdf")

# ---- Hallmark IFN-α Response Over Time ----

library(GSVA)
library(msigdbr)
library(ggplot2)
library(dplyr)

# Load saved objects
v <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/voom(hSpS_Ctrx_batch).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx_batch).rds")

# Extract log-CPM normalized expression matrix
norm_mat <- v$E

# Retrieve Hallmark Interferon Alpha Response gene set
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
ifna_genes <- hallmark %>%
  filter(gs_name == "HALLMARK_INTERFERON_ALPHA_RESPONSE") %>%
  pull(gene_symbol) %>%
  unique()

gene_sets <- list(HALLMARK_INTERFERON_ALPHA_RESPONSE = ifna_genes)

# Compute per-sample ssGSEA enrichment scores
gsva_params <- ssgseaParam(exprData = norm_mat, geneSets = gene_sets)
scores <- gsva(gsva_params)

# Combine scores with sample metadata
score_df <- data.frame(
  Admera_Health_ID = colnames(scores),
  ifna_score = as.numeric(scores["HALLMARK_INTERFERON_ALPHA_RESPONSE", ])
) %>%
  left_join(meta_filtered, by = "Admera_Health_ID")

# Summarize mean and standard error per group and timepoint
summary_df <- score_df %>%
  group_by(Line, Stage) %>%
  summarise(
    mean_score = mean(ifna_score),
    se = sd(ifna_score) / sqrt(n()),
    .groups = "drop"
  )

# Line plot of pathway score over organoid age
ggplot(summary_df, aes(x = Stage, y = mean_score, color = Line, group = Line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.15) +
  labs(
    title = "Hallmark Interferon Alpha Response Over Time",
    x = "Organoid Age",
    y = "ssGSEA Enrichment Score",
    color = "Group"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_manual(values = c("control" = "#2171B5", "ALS" = "#E66101"))

# ---- Hallmark Pathway Scores Over Time (Two-Panel) ----

library(GSVA)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load saved objects
v <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/voom(hSpS_Ctrx_batch).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx_batch).rds")

# Extract log-CPM normalized expression matrix
norm_mat <- v$E

# Define pathways of interest
pathway_names <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

# Retrieve gene sets from msigdbr
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- hallmark %>%
  filter(gs_name %in% pathway_names) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

# Compute per-sample ssGSEA enrichment scores
gsva_params <- ssgseaParam(exprData = norm_mat, geneSets = gene_sets)
scores <- gsva(gsva_params)

# Reshape scores to long format and join metadata
score_df <- as.data.frame(t(scores)) %>%
  mutate(Admera_Health_ID = rownames(.)) %>%
  pivot_longer(-Admera_Health_ID, names_to = "pathway", values_to = "score") %>%
  left_join(meta_filtered, by = "Admera_Health_ID")

# Clean pathway names for display
score_df$pathway <- score_df$pathway %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

# Summarize mean and standard error per group, timepoint, and pathway
summary_df <- score_df %>%
  group_by(Line, Stage, pathway) %>%
  summarise(
    mean_score = mean(score),
    se = sd(score) / sqrt(n()),
    .groups = "drop"
  )

# Faceted by pathway, control vs ALS lines in each panel
ggplot(summary_df, aes(x = Stage, y = mean_score, color = Line, group = Line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.15) +
  facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  labs(
    x = "Organoid Age",
    y = "ssGSEA Enrichment Score",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("control" = "#2171B5", "ALS" = "#E66101"))


# ---- Hallmark Pathway Scores Over Time ----
library(GSVA)
library(msigdbr)
library(limma)
library(ggplot2)
library(dplyr)
library(tidyr)

# Load saved objects
v <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/voom(hSpS_Ctrx_batch).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx_batch).rds")
fit2 <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/fit2(hSpS_Ctrx_batch).rds")

gsea_d30  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d30_ALS_vs_Ctrl(Batch).rds")
gsea_d50  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d50_ALS_vs_Ctrl(Batch).rds")
gsea_d75  <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d75_ALS_vs_Ctrl(Batch).rds")
gsea_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d120_ALS_vs_Ctrl(Batch).rds")

# Remove batch and sex effects from expression matrix
norm_mat <- removeBatchEffect(v$E,
                              batch = meta_filtered$Batch,
                              covariates = model.matrix(~ Sex, data = meta_filtered)[, -1]
)

# Define pathways of interest
pathway_names <- c(
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

# Retrieve gene sets from msigdbr
hallmark <- msigdbr(species = "Homo sapiens", category = "H")
gene_sets <- hallmark %>%
  filter(gs_name %in% pathway_names) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

# Compute per-sample ssGSEA enrichment scores
gsva_params <- ssgseaParam(exprData = norm_mat, geneSets = gene_sets)
scores <- gsva(gsva_params)

# Reshape scores to long format and join metadata
score_df <- as.data.frame(t(scores)) %>%
  mutate(Admera_Health_ID = rownames(.)) %>%
  pivot_longer(-Admera_Health_ID, names_to = "pathway", values_to = "score") %>%
  left_join(meta_filtered, by = "Admera_Health_ID")

# Clean pathway names for display
score_df$pathway <- score_df$pathway %>%
  gsub("HALLMARK_", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

# Summarize mean and standard error per group, timepoint, and pathway
summary_df <- score_df %>%
  group_by(Line, Stage, pathway) %>%
  summarise(
    mean_score = mean(score),
    se = sd(score) / sqrt(n()),
    .groups = "drop"
  )

# Combine pre-computed fgsea results and tag with stage
fgsea_results <- bind_rows(
  gsea_d30  %>% mutate(Stage = "d30"),
  gsea_d50  %>% mutate(Stage = "d50"),
  gsea_d75  %>% mutate(Stage = "d75"),
  gsea_d120 %>% mutate(Stage = "d120")
)

# Filter to pathways of interest and add significance labels
stat_df <- fgsea_results %>%
  filter(pathway %in% pathway_names) %>%
  mutate(
    pathway = pathway %>%
      gsub("HALLMARK_", "", .) %>%
      gsub("_", " ", .) %>%
      tools::toTitleCase(),
    Stage = factor(Stage, levels = c("d30", "d50", "d75", "d120")),
    label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  select(pathway, Stage, padj, label)

# Position stars above the higher group at each timepoint
star_pos <- summary_df %>%
  group_by(pathway, Stage) %>%
  summarise(y_pos = max(mean_score + se) + 0.02, .groups = "drop") %>%
  left_join(stat_df, by = c("pathway", "Stage"))

# Faceted line plot with significance stars
ggplot(summary_df, aes(x = Stage, y = mean_score, color = Line, group = Line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.15) +
  geom_text(data = star_pos, aes(x = Stage, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  labs(
    x = "Organoid Age",
    y = "ssGSEA Enrichment Score",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("control" = "#2171B5", "ALS" = "#E66101"))

ggsave(
  filename = "figures/Shilu's vs Jimena/hallmark_pathway_scores_over_time.pdf",
  width = 16,
  height = 5,
  dpi = 300
)




# ── ISG Heatmap ──────────────────────────────────────────────────────────────
library(limma)
library(pheatmap)
library(dplyr)

# Define core ISG panel
isg_genes <- c(
  "STAT1", "STAT2", "IRF7", "IRF9",
  "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
  "ISG15", "ISG20", "RSAD2", "IFIT1", "IFIT2", "IFIT3",
  "IFITM1", "IFITM3", "CXCL10", "CXCL11", "CCL5",
  "BST2", "IFI44", "IFI44L", "IFI6", "IFI27",
  "HERC5", "USP18", "DDX58", "IFIH1", "TRIM22"
)

# Batch-correct log2-CPM
design_interest <- model.matrix(~ Line, data = meta_filtered)
norm_mat <- removeBatchEffect(v$E,
                              batch = meta_filtered$Batch,
                              covariates = model.matrix(~ Sex, data = meta_filtered)[, -1],
                              design = design_interest)

# Subset to ISGs present in the data
isg_present <- intersect(isg_genes, rownames(norm_mat))
mat <- norm_mat[isg_present, ]

# Z-score per timepoint
stages <- c("d30", "d50", "d75", "d120")
mat_z_list <- list()
for (stg in stages) {
  idx <- meta_filtered$Admera_Health_ID[meta_filtered$Stage == stg]
  idx <- intersect(idx, colnames(mat))
  sub <- mat[, idx]
  mat_z_list[[stg]] <- t(scale(t(sub)))
}
mat_z <- do.call(cbind, mat_z_list)

# Order: Stage -> Line -> cluster within each Stage x Line block
groups <- c("control", "ALS")
ordered_ids <- c()
for (stg in stages) {
  for (grp in groups) {
    idx <- meta_filtered$Admera_Health_ID[meta_filtered$Stage == stg & meta_filtered$Line == grp]
    idx <- intersect(idx, colnames(mat_z))
    if (length(idx) > 2) {
      hc <- hclust(dist(t(mat_z[, idx])))
      idx <- idx[hc$order]
    }
    ordered_ids <- c(ordered_ids, idx)
  }
}
mat_z <- mat_z[, ordered_ids]

# Column annotation
col_ann <- meta_filtered %>%
  filter(Admera_Health_ID %in% colnames(mat_z)) %>%
  as.data.frame()
rownames(col_ann) <- col_ann$Admera_Health_ID
col_ann <- col_ann[ordered_ids, c("Line", "Stage")]

# Gaps between timepoints
gap_positions <- cumsum(table(col_ann$Stage))
gap_positions <- gap_positions[-length(gap_positions)]

# Annotation colors
ann_colors <- list(
  Line  = c("control" = "#2171B5", "ALS" = "#E66101"),
  Stage = c("d30" = "#BDBDBD", "d50" = "#969696", "d75" = "#636363", "d120" = "#252525")
)

# Heatmap
pheatmap(mat_z,
         annotation_col = col_ann,
         annotation_colors = ann_colors,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_colnames = FALSE,
         scale = "none",
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(-2, 2, length.out = 101),
         gaps_col = as.numeric(gap_positions),
         border_color = NA,
         fontsize_row = 8,
         main = "ISG Expression — ALS vs Control")

# ── ISG / Neuronal Heatmap ───────────────────────────────────────────────────
library(limma)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# Gene panels
isg_genes <- c(
  "STAT1", "STAT2", "IRF7", "IRF9",
  "MX1", "MX2", "OAS1", "OAS2", "OAS3", "OASL",
  "ISG15", "ISG20", "RSAD2", "IFIT1", "IFIT2", "IFIT3",
  "IFITM1", "IFITM3", "CXCL10", "CXCL11", "CCL5",
  "BST2", "IFI44", "IFI44L", "IFI6", "IFI27",
  "HERC5", "USP18", "DDX58", "IFIH1", "TRIM22"
)
neuronal_genes     <- c("CAMK2A", "STMN2", "RBFOX3", "NEFL", "THY1", "MAPT", "MAP2", "TUBB3")
synaptic_genes     <- c("SNAP25", "SYN1", "SYN2", "SYT1", "SYT2", "SYP")
motor_neuron_genes <- c("ISL1", "SLC5A7", "SLC18A3", "CHAT", "PHOX2B", "MNX1")

all_genes <- unique(c(isg_genes, neuronal_genes, synaptic_genes, motor_neuron_genes))

# Batch-correct log2-CPM
design_interest <- model.matrix(~ Line, data = meta_filtered)
norm_mat <- removeBatchEffect(v$E,
                              batch = meta_filtered$Batch,
                              covariates = model.matrix(~ Sex, data = meta_filtered)[, -1],
                              design = design_interest)

# Subset to genes present in the data
genes_present <- intersect(all_genes, rownames(norm_mat))
mat <- norm_mat[genes_present, ]

# Row z-score
mat_z <- t(scale(t(mat)))

# Order columns: Stage -> Line -> cluster within each block
stages <- c("d30", "d50", "d75", "d120")
groups <- c("control", "ALS")
ordered_ids <- c()
for (stg in stages) {
  for (grp in groups) {
    idx <- meta_filtered$Admera_Health_ID[meta_filtered$Stage == stg & meta_filtered$Line == grp]
    idx <- intersect(idx, colnames(mat_z))
    if (length(idx) > 2) {
      hc <- hclust(dist(t(mat_z[, idx])))
      idx <- idx[hc$order]
    }
    ordered_ids <- c(ordered_ids, idx)
  }
}
mat_z <- mat_z[, ordered_ids]

# Column split factor
col_meta <- meta_filtered %>%
  filter(Admera_Health_ID %in% ordered_ids) %>%
  as.data.frame()
rownames(col_meta) <- col_meta$Admera_Health_ID
col_meta <- col_meta[ordered_ids, ]
col_split <- factor(paste0(col_meta$Stage, " - ", col_meta$Line),
                    levels = paste0(rep(stages, each = 2), " - ", rep(groups, times = 4)))

# Row annotation
row_cat <- case_when(
  genes_present %in% isg_genes          ~ "ISG",
  genes_present %in% neuronal_genes     ~ "General Neuronal",
  genes_present %in% synaptic_genes     ~ "Synaptic",
  genes_present %in% motor_neuron_genes ~ "Motor Neuron"
)

row_ha <- rowAnnotation(
  Category = row_cat,
  col = list(Category = c("ISG" = "#7570B3", "General Neuronal" = "#1B9E77",
                          "Synaptic" = "#D95F02", "Motor Neuron" = "#E7298A")),
  annotation_legend_param = list(
    Category = list(labels_rot = 0)
  )
)

# Column annotations
col_ha <- HeatmapAnnotation(
  Line  = col_meta$Line,
  Stage = col_meta$Stage,
  col = list(
    Line  = c("control" = "#2171B5", "ALS" = "#E66101"),
    Stage = c("d30" = "#BDBDBD", "d50" = "#969696", "d75" = "#636363", "d120" = "#252525")
  )
)

# Color scale
col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B"))

# Build heatmap object
ht <- Heatmap(mat_z,
              name = "Z-score",
              col = col_fun,
              top_annotation = col_ha,
              right_annotation = row_ha,
              column_split = col_split,
              column_gap = unit(c(1, 3, 1, 3, 1, 3, 1), "mm"),
              column_title = NULL,
              cluster_columns = TRUE,
              cluster_column_slices = FALSE,
              cluster_rows = TRUE,
              show_column_names = FALSE,
              show_row_names = TRUE,
              row_names_gp = gpar(fontsize = 6),
              border = TRUE,
              heatmap_legend_param = list(title = "Z-score"))

# Save to PDF
pdf("figures/Shilu's vs Jimena/isg_neuronal_heatmap.pdf", width = 14, height = 10)
draw(ht, column_title = "ISG & Neuronal Marker Expression — ALS vs Control",
     column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

# ── GOBP Neuronal ssGSEA with fgsea stats ───────────────────────────────────
library(limma)
library(GSVA)
library(msigdbr)
library(fgsea)
library(dplyr)
library(tidyr)
library(ggplot2)

# Batch-correct log2-CPM
design_interest <- model.matrix(~ Line, data = meta_filtered)
norm_mat <- removeBatchEffect(v$E,
                              batch = meta_filtered$Batch,
                              covariates = model.matrix(~ Sex, data = meta_filtered)[, -1],
                              design = design_interest)

# Define neuronal pathways
pathway_names <- c(
  "GOBP_SYNAPTIC_SIGNALING",
  "GOBP_NEUROTRANSMITTER_TRANSPORT",
  "GOBP_MOTOR_NEURON_AXON_GUIDANCE",
  "GOBP_NEURON_APOPTOTIC_PROCESS",
  "GOBP_SYNAPSE_ORGANIZATION"
)

# Retrieve gene sets
gobp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")
gene_sets <- gobp %>%
  filter(gs_name %in% pathway_names) %>%
  split(.$gs_name) %>%
  lapply(function(x) unique(x$gene_symbol))

# Compute ssGSEA enrichment scores
gsva_params <- ssgseaParam(exprData = norm_mat, geneSets = gene_sets)
scores <- gsva(gsva_params)

# Reshape to long format and join metadata
score_df <- as.data.frame(t(scores)) %>%
  mutate(Admera_Health_ID = rownames(.)) %>%
  pivot_longer(-Admera_Health_ID, names_to = "pathway", values_to = "score") %>%
  left_join(meta_filtered, by = "Admera_Health_ID")

# Clean pathway names for display
score_df$pathway <- score_df$pathway %>%
  gsub("GOBP_", "", .) %>%
  gsub("_", " ", .) %>%
  tools::toTitleCase()

# Summarize mean and SE per group, timepoint, and pathway
summary_df <- score_df %>%
  group_by(Line, Stage, pathway) %>%
  summarise(
    mean_score = mean(score),
    se = sd(score) / sqrt(n()),
    .groups = "drop"
  )

# Run fgsea using moderated t-stats from eBayes
contrasts <- c("d30_ALS", "d50_ALS", "d75_ALS", "d120_ALS")
stage_labels <- c("d30", "d50", "d75", "d120")

fgsea_results <- lapply(seq_along(contrasts), function(i) {
  tt <- topTable(fit2, coef = contrasts[i], number = Inf, sort.by = "none")
  stats <- setNames(tt$t, rownames(tt))
  stats <- sort(stats, decreasing = TRUE)
  
  res <- fgsea(pathways = gene_sets, stats = stats, minSize = 15, maxSize = 500, nproc = 4)
  res$Stage <- stage_labels[i]
  res
}) %>%
  bind_rows()

# Filter and add significance labels
stat_df <- fgsea_results %>%
  mutate(
    pathway = pathway %>%
      gsub("GOBP_", "", .) %>%
      gsub("_", " ", .) %>%
      tools::toTitleCase(),
    Stage = factor(Stage, levels = c("d30", "d50", "d75", "d120")),
    label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ ""
    )
  ) %>%
  select(pathway, Stage, padj, label)

# Position stars above the higher group
star_pos <- summary_df %>%
  group_by(pathway, Stage) %>%
  summarise(y_pos = max(mean_score + se) + 0.02, .groups = "drop") %>%
  left_join(stat_df, by = c("pathway", "Stage"))

# Faceted line plot with significance stars
ggplot(summary_df, aes(x = Stage, y = mean_score, color = Line, group = Line)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se), width = 0.15) +
  geom_text(data = star_pos, aes(x = Stage, y = y_pos, label = label),
            inherit.aes = FALSE, size = 5, color = "black") +
  facet_wrap(~ pathway, scales = "free_y", nrow = 1) +
  labs(
    x = "Organoid Age",
    y = "ssGSEA Enrichment Score",
    color = "Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_color_manual(values = c("control" = "#2171B5", "ALS" = "#E66101"))

ggsave(
  filename = "figures/Shilu's vs Jimena/gobp_neuronal_ssgsea_over_time.pdf",
  width = 16,
  height = 5,
  dpi = 300
)
# ---- D120 gsea and volcano ----
library(tidyverse)
library(ggrepel)
library(patchwork)

gsea_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d120_ALS_vs_Ctrl(Batch).rds")
# Load your per-stage limma results
res_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/res_d120_ALS_vs_Ctrl(Batch).rds")

# Panel 1: Lollipop bubble plot of significant hallmark pathways ---

sig_paths <- gsea_d120 %>%
  filter(abs(NES) > 1.5, pval < 0.05) %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    direction = ifelse(NES > 0, "Enriched", "Suppressed"),
    neg_log10_padj = -log10(padj)
  ) %>%
  arrange(NES)

sig_paths$pathway_clean <- factor(sig_paths$pathway_clean, levels = sig_paths$pathway_clean)

p1 <- ggplot(sig_paths, aes(x = NES, y = pathway_clean)) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean),
               color = "grey60", linewidth = 0.4) +
  geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("Enriched" = "#D73027", "Suppressed" = "#4575B4")) +
  scale_size_continuous(range = c(2, 7), name = "-log10(p-value)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Direction",
       title = "C9orf72 vs. Control (Day 120)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))
p1

# Extract leading edge genes from innate immunity-related pathway(s)
innate_paths <- gsea_d120 %>%
  filter(str_detect(pathway, regex("INNATE_IMMUN|INFLAMMATORY_RESPONSE|INTERFERON", ignore_case = TRUE)))
innate_genes <- unique(unlist(innate_paths$leadingEdge))

volcano_df <- res_d120 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    FC = 2^logFC,
    log10_p = log10(P.Value),
    is_innate = gene %in% innate_genes,
    innate_direction = case_when(
      is_innate & logFC > 0 ~ "Up",
      is_innate & logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  )

innate_min_p <- volcano_df %>%
  filter(is_innate) %>%
  pull(log10_p) %>%
  min(na.rm = TRUE)

y_floor <- innate_min_p * 1.1


p2 <- ggplot(volcano_df, aes(x = FC, y = log10_p)) +
  geom_point(data = filter(volcano_df, !is_innate),
             color = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = filter(volcano_df, is_innate),
             aes(color = innate_direction), size = 1.8, alpha = 0.9) +
  geom_text_repel(data = filter(volcano_df, is_innate),
                  aes(label = gene), size = 2.8, fontface = "italic",
                  max.overlaps = 30, segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.3,
                  xlim = c(1/4, 4), ylim = c(y_floor, 0)) +
  scale_color_manual(values = c("Up" = "red", "Down" = "black"), guide = "none") +
  scale_x_continuous(
    trans = "log2",
    breaks = c(1/4, 1/2, 1, 2, 4),
    labels = c("-4", "-2", "1", "+2", "+4")
  ) +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(x = "FC: C9orf72 vs. Control", y = "-log10(p-value)",
       title = "Innate Immunity Genes (Day 120)") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  ) +
  coord_cartesian(xlim = c(1/4, 4), ylim = c(y_floor, 0), clip = "on")
p2

# --- Combine ---
combined <- p1 + p2 +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A")

combined
ggsave("figures/Shilu's vs Jimena/d120_volcano.pdf", p2, width = 7, height = 5, dpi = 300)

# ---- Generalizable GSEA and Volcano Jimena----
# -- User parameters
day_label     <- "d120"
comparison    <- "ALS_vs_Ctrl(Batch)"
pathway_regex <- "OXIDATIVE_PHOSPHORYLATION"
nes_cutoff    <- 1.5
pval_cutoff   <- 0.05
b_title       <- "Oxidative Phosphorylation Genes"

# -- Load data 
base_dir     <- "data/Anderson-Suthar_collab/r_objects"
gsea_results <- readRDS(file.path(base_dir, paste0("fgsea_", day_label, "_", comparison, ".rds")))
res_results  <- readRDS(file.path(base_dir, paste0("res_", day_label, "_", comparison, ".rds")))

day_display <- str_replace(day_label, "^d", "Day ")

# Filter to significant pathways and clean names for display
sig_paths <- gsea_results %>%
  filter(abs(NES) > nes_cutoff, pval < pval_cutoff) %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_") %>%
      str_replace_all("_", " ") %>%
      str_to_title(),
    direction = ifelse(NES > 0, "Enriched", "Suppressed"),
    neg_log10_padj = -log10(padj)
  ) %>%
  arrange(NES)

sig_paths$pathway_clean <- factor(sig_paths$pathway_clean, levels = sig_paths$pathway_clean)

# Lollipop bubble plot of significant hallmark pathways
p1 <- ggplot(sig_paths, aes(x = NES, y = pathway_clean)) +
  geom_segment(aes(x = 0, xend = NES, y = pathway_clean, yend = pathway_clean),
               color = "grey60", linewidth = 0.4) +
  geom_point(aes(size = neg_log10_padj, fill = direction), shape = 21, stroke = 0.3) +
  scale_fill_manual(values = c("Enriched" = "#D73027", "Suppressed" = "#4575B4")) +
  scale_size_continuous(range = c(2, 7), name = "-log10(p-value)") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL, fill = "Direction",
       title = paste0("C9orf72 vs. Control (", day_display, ")")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major.y = element_blank(), legend.position = "right",
        plot.title = element_text(face = "bold", size = 12))

# Pull leading-edge genes from pathways matching the target regex
target_paths <- gsea_results %>%
  filter(str_detect(pathway, regex(pathway_regex, ignore_case = TRUE)))
target_genes <- unique(unlist(target_paths$leadingEdge))

# Build volcano data frame with fold change and target gene flags
volcano_df <- res_results %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    FC = 2^logFC,
    log10_p = log10(P.Value),
    is_target = gene %in% target_genes,
    target_direction = case_when(
      is_target & logFC > 0 ~ "Up",
      is_target & logFC < 0 ~ "Down",
      TRUE ~ NA_character_
    )
  )

# Set y-axis floor based on the least significant target gene
target_min_p <- volcano_df %>%
  filter(is_target) %>%
  pull(log10_p) %>%
  min(na.rm = TRUE)
y_floor <- target_min_p * 1.2

# Readable title derived from the pathway regex
pathway_title <- str_replace_all(pathway_regex, "\\|", " / ") %>%
  str_replace_all("_", " ") %>%
  str_to_title()

# Volcano plot highlighting leading-edge genes
p2 <- ggplot(volcano_df, aes(x = FC, y = log10_p)) +
  geom_point(data = filter(volcano_df, !is_target),
             color = "grey70", size = 0.6, alpha = 0.4) +
  geom_point(data = filter(volcano_df, is_target),
             aes(color = target_direction), size = 1.8, alpha = 0.9) +
  geom_text_repel(data = filter(volcano_df, is_target),
                  aes(label = gene), size = 2.8, fontface = "italic",
                  max.overlaps = 30, segment.size = 0.3,
                  min.segment.length = 0, box.padding = 0.3) +
  scale_color_manual(values = c("Up" = "red", "Down" = "black"), guide = "none") +
  scale_x_continuous(
    trans = "log2",
    breaks = c(1/4, 1/2, 1, 2, 4),
    labels = c("-4", "-2", "1", "+2", "+4")
  ) +
  coord_cartesian(xlim = c(1/4, 4), ylim = c(y_floor, 0), clip = "on") +
  geom_vline(xintercept = 1, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  labs(x = "FC: C9orf72 vs. Control", y = "-log10(p-value)",
       title = paste0(b_title, " (", day_display, ")")) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

# Combine panels side by side
combined <- p1 + p2 +
  plot_layout(widths = c(1, 1.5)) +
  plot_annotation(tag_levels = "A")

combined
ggsave("figures/Shilu's vs Jimena/d120_volcano.pdf", p2, width = 7, height = 5, dpi = 300)
####

# ---- Innate Heatmap ----
library(tidyverse)
library(pheatmap)

gsea_d120 <- readRDS("data/Anderson-Suthar_collab/r_objects/fgsea_d120_ALS_vs_Ctrl(Batch).rds")
v <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/voom(hSpS_Ctrx_batch).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx_batch).rds")

# Get innate immunity leading edge genes
innate_paths <- gsea_d120 %>%
  filter(str_detect(pathway, regex("INNATE_IMMUN|INFLAMMATORY_RESPONSE|INTERFERON", ignore_case = TRUE)))
innate_genes <- unique(unlist(innate_paths$leadingEdge))

# Subset to d120
d120_idx <- meta_filtered$Stage == "d120"
expr_d120 <- v$E[, d120_idx]
meta_d120 <- meta_filtered[d120_idx, ]

# Check actual Line levels and subset
print(levels(meta_d120$Line))
print(table(meta_d120$Line))

# Subset to innate genes
innate_expr <- expr_d120[rownames(expr_d120) %in% innate_genes, ]

# Center relative to control — using case-insensitive matching
ctrl_idx <- grepl("control", meta_d120$Line, ignore.case = TRUE)
ctrl_means <- rowMeans(innate_expr[, ctrl_idx])
innate_centered <- innate_expr - ctrl_means

# Annotation
anno_col <- data.frame(
  Line = meta_d120$Line,
  row.names = colnames(innate_centered)
)

# Build color mapping dynamically from actual levels
line_levels <- levels(meta_d120$Line)
line_colors <- setNames(
  c("#D73027", "#4575B4")[seq_along(line_levels)],
  line_levels
)
# Order columns: control first, then ALS
col_order <- order(meta_d120$Line)
innate_ordered <- innate_centered[, col_order]
anno_col_ordered <- data.frame(
  Line = meta_d120$Line[col_order],
  row.names = colnames(innate_ordered)
)

pheatmap(innate_ordered,
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         show_colnames = FALSE,
         annotation_col = anno_col_ordered,
         annotation_colors = anno_colors,
         color = colorRampPalette(c("#4575B4", "white", "#D73027"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color = NA,
         fontsize_row = 5,
         main = "Innate Immunity Genes Relative to Control (Day 120)")
