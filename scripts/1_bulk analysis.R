# ---- Library ----

library(tidyverse)
library(DESeq2)
library(rtracklayer)
library(apeglm)
library(pheatmap)

# ---- Data import ----
# GTF for gene_ids
gtf <- import("data/gtf/Homo_sapiens.GRCh38.114.gtf")
gtf_genes <- gtf[gtf$type == "gene"]

# Import all count files
count_files <- list.files(
  "data/rnaseq/trimmed_aligned",
  pattern = "_featureCounts_exon.txt$",
  recursive = TRUE,
  full.names = TRUE
)
count_files

counts_list <- lapply(count_files, function(f) {
  df <- read.delim(f, comment.char = "#", check.names = FALSE)
  
  counts <- df[, c(1, ncol(df))]
  sample <- sub("_featureCounts_exon.txt$", "", basename(f))
  colnames(counts) <- c("gene_id", sample)
  
  counts
})

#Merge count files into one matrix
counts_merged <- Reduce(function(x, y) merge(x, y, by = "gene_id", all = TRUE), counts_list)
counts_merged[is.na(counts_merged)] <- 0

rownames(counts_merged) <- counts_merged$gene_id
counts_merged <- counts_merged[, -1]


#Add gene_id column for reference
gene_map <- setNames(
  gtf_genes$gene_name,
  gtf_genes$gene_id
)

counts_merged$gene_name <- gene_map[rownames(counts_merged)]

write.csv(counts_merged, "data/R2SDHF_counts_per_sample(gene_id).csv", quote = FALSE)

# ---- Build meta Data ----
#create seperate matrix without gene_id
counts_mat <- counts_merged |>
  dplyr::select(-gene_name) |>
  as.matrix()

mode(counts_mat) <- "integer"

#Create sample names and meta data
sample_names <- colnames(counts_mat)

colData <- data.frame(
  sample = sample_names
)

colData <- colData |>
  mutate(
    group = case_when(
      str_detect(sample, "Mock") ~ "Mock",
      TRUE ~ str_extract(sample, "TX|MAD|ROV")
    ),
    time = case_when(
      str_detect(sample, "Mock") ~ "0",
      TRUE ~ str_extract(sample, "12|24|48")
    ),
    condition = ifelse(group == "Mock",
                       "Mock",
                       paste0(group, "_", time))
  )

rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, "Mock")


# ---- Run DESeq2 ----
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat[rowSums(counts_mat) > 0, ],
  colData   = colData,
  design    = ~ condition
)

dds <- DESeq(dds)

# ---- Extract comps relative to Mock ----
res_list <- lapply(
  levels(colData$condition)[levels(colData$condition) != "Mock"],
  function(cond) {
    res <- results(dds, contrast = c("condition", cond, "Mock"))
    res_df <- as.data.frame(res)
    res_df$ensembl_id <- rownames(res_df)
    res_df
  }
)

names(res_list) <- levels(colData$condition)[levels(colData$condition) != "Mock"]

#compare results without shrinkage
res <- results(dds)
plotMA(res)
plotMA(res)

# ---- Shrink data ----
resultsNames(dds)

coef_names <- resultsNames(dds)
coef_names <- coef_names[grep("vs_Mock", coef_names)]

res_shrunk_list <- lapply(coef_names, function(coef_name) {
  res <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  res_df
})

names(res_shrunk_list) <- coef_names

saveRDS(res_shrunk_list, "data/R2SDHF/res_shrunk_list.rds")

# ---- PCA with all samples-----
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data$treatment_time <- as.character(pca_data$condition)

base_cols <- c(
  Mock = "#D62728",
  TX   = "blue",
  MAD  = "darkgreen",
  ROV  = "purple"
)

make_shades <- function(hex, n = 3) {
  colorRampPalette(c("#FFFFFF", hex))(n + 1)[-1]
}

tx_shades  <- make_shades(base_cols[["TX"]],  3)
mad_shades <- make_shades(base_cols[["MAD"]], 3)
rov_shades <- make_shades(base_cols[["ROV"]], 3)

col_map <- c(
  Mock   = base_cols[["Mock"]],
  TX_12  = tx_shades[1],  TX_24  = tx_shades[2],  TX_48  = tx_shades[3],
  MAD_12 = mad_shades[1], MAD_24 = mad_shades[2], MAD_48 = mad_shades[3],
  ROV_12 = rov_shades[1], ROV_24 = rov_shades[2], ROV_48 = rov_shades[3]
)

pca_data$treatment_time <- factor(
  pca_data$treatment_time,
  levels = c("Mock","TX_12","TX_24","TX_48","MAD_12","MAD_24","MAD_48","ROV_12","ROV_24","ROV_48")
)

centroids <- pca_data %>%
  group_by(treatment_time) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

p <- ggplot(pca_data, aes(PC1, PC2, color = treatment_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(
    data = centroids,
    aes(PC1, PC2, fill = treatment_time),
    shape = 24,           # filled triangle
    color = "black",
    size = 3.5,
    stroke = 1.2
  ) +
  scale_color_manual(values = col_map, drop = FALSE) +
  scale_fill_manual(values = col_map, drop = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()

ggsave(
  filename = "figures/PCA_treatment_time.png",
  plot = p,
  width = 8,
  height = 6
)

# ---- Volcano plots ----
padj_cutoff <- 0.05
lfc_cutoff  <- 2

outdir <- "figures/volcano"

for (coef_name in names(res_shrunk_list)) {
  
  volcano_df <- res_shrunk_list[[coef_name]] |>
    mutate(
      neglog10_padj = -log10(padj),
      direction = case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  n_up   <- sum(volcano_df$direction == "Up", na.rm = TRUE)
  n_down <- sum(volcano_df$direction == "Down", na.rm = TRUE)
  
  x_max <- max(volcano_df$log2FoldChange, na.rm = TRUE)
  x_min <- min(volcano_df$log2FoldChange, na.rm = TRUE)
  y_max <- max(volcano_df$neglog10_padj, na.rm = TRUE)
  
  p <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10_padj)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    annotate("text",
             x = x_max, y = y_max,
             label = paste0("Up: ", n_up),
             hjust = 1, vjust = 1,
             color = "red") +
    annotate("text",
             x = x_min, y = y_max,
             label = paste0("Down: ", n_down),
             hjust = 0, vjust = 1,
             color = "blue") +
    scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    theme_classic() +
    labs(
      title = coef_name,
      x = "Shrunken log2 Fold Change",
      y = "-log10 adjusted p-value"
    )
  
  comp_label <- sub("^condition_", "", coef_name)
  comp_label <- sub("_vs_Mock$", "_vs_Mock", comp_label)
  
  outfile <- file.path(outdir, paste0("volcano_", comp_label, ".png"))
  
  ggsave(outfile, p, width = 7, height = 6, dpi = 300)
}

#### One Example ####
# padj_cutoff <- 0.05
# lfc_cutoff  <- 2
# 
# coef_name <- names(res_shrunk_list)[1]
# 
# volcano_df <- res_shrunk_list[[coef_name]] |>
#   mutate(
#     neglog10_padj = -log10(padj),
#     direction = case_when(
#       padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
#       padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
#       TRUE ~ "NS"
#     )
#   )
# 
# n_up   <- sum(volcano_df$direction == "Up", na.rm = TRUE)
# n_down <- sum(volcano_df$direction == "Down", na.rm = TRUE)
# 
# x_max <- max(volcano_df$log2FoldChange, na.rm = TRUE)
# x_min <- min(volcano_df$log2FoldChange, na.rm = TRUE)
# y_max <- max(volcano_df$neglog10_padj, na.rm = TRUE)
# 
# p <- ggplot(volcano_df, aes(x = log2FoldChange, y = neglog10_padj)) +
#   geom_point(aes(color = direction), alpha = 0.6, size = 1) +
#   geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
#   geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
#   annotate("text",
#            x = x_max, y = y_max,
#            label = paste0("Up: ", n_up),
#            hjust = 1, vjust = 1,
#            color = "red") +
#   annotate("text",
#            x = x_min, y = y_max,
#            label = paste0("Down: ", n_down),
#            hjust = 0, vjust = 1,
#            color = "blue") +
#   scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
#   theme_classic() +
#   labs(
#     title = coef_name,
#     x = "Shrunken log2 Fold Change",
#     y = "-log10 adjusted p-value"
#   )
# 
# p
# 
# 

# ---- Heatmap ----
# VST
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

# Ensembl -> gene_name map (data.frame, not a vector)
gene_map <- data.frame(
  ensembl_id = rownames(counts_merged),
  gene_name  = counts_merged$gene_name,
  stringsAsFactors = FALSE
)

# Match map to VST rows
idx <- match(rownames(mat), gene_map$ensembl_id)
mapped_names <- gene_map$gene_name[idx]

# Top 50 most variable genes
gene_vars <- apply(mat, 1, var)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:50]

mat_subset <- mat[top_genes, ]

# Replace Ensembl rownames with gene_name (fallback to Ensembl if missing/blank)
new_names <- mapped_names[match(top_genes, rownames(mat))]
new_names[is.na(new_names) | new_names == ""] <- top_genes[is.na(new_names) | new_names == ""]
rownames(mat_subset) <- new_names

# Z-score per gene
mat_scaled <- t(scale(t(mat_subset)))

# Parse sample names: R2SDHF_<id>_<treatment>_<time>_...
sample_names <- colnames(mat_scaled)
parts <- strsplit(sample_names, "_", fixed = TRUE)

treatment <- vapply(parts, function(x) x[3], character(1))
time <- vapply(parts, function(x) x[4], character(1))

treatment[treatment == "Mock"] <- "Mock"
time[treatment == "Mock"] <- "0"

annotation_col <- data.frame(
  treatment = treatment,
  time = time,
  row.names = sample_names,
  stringsAsFactors = FALSE
)

# Order columns: Mock first, then TX, MAD, ROV by time
desired_order <- annotation_col |>
  mutate(
    treatment = factor(treatment, levels = c("Mock","TX","MAD","ROV")),
    time = as.numeric(time)
  ) |>
  arrange(treatment, time) |>
  rownames()

mat_scaled <- mat_scaled[, desired_order]
annotation_col <- annotation_col[desired_order, ]

# Heatmap
pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE
)