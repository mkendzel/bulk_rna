# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)

# ----- Import data ----
expr <- read.delim(
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- set up data frame ----
# Identify count columns
count_cols <- grep("_count$", colnames(expr), value = TRUE)

# Subset
counts <- expr[, c("gene_id", count_cols)]

# Set rownames
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# Remove "_count" suffix
colnames(counts) <- sub("_count$", "", colnames(counts))


# rename columns
# Map old -> new (explicit)
rename_map <- c(
  "R2SDHF_1"  = "Mock_0_1",
  "R2SDHF_2"  = "Mock_0_2",
  "R2SDHF_3"  = "Mock_0_3",
  "R2SDHF_4"  = "TX_12_1",
  "R2SDHF_5"  = "TX_12_2",
  "R2SDHF_6"  = "TX_12_3",
  "R2SDHF_7"  = "ROV_12_1",
  "R2SDHF_8"  = "ROV_12_2",
  "R2SDHF_9"  = "ROV_12_3",
  "R2SDHF_10" = "MAD_12_1",
  "R2SDHF_11" = "MAD_12_2",
  "R2SDHF_12" = "MAD_12_3",
  "R2SDHF_13" = "TX_24_1",
  "R2SDHF_14" = "TX_24_2",
  "R2SDHF_15" = "TX_24_3",
  "R2SDHF_16" = "ROV_24_1",
  "R2SDHF_17" = "ROV_24_2",
  "R2SDHF_18" = "ROV_24_3",
  "R2SDHF_19" = "MAD_24_1",
  "R2SDHF_20" = "MAD_24_2",
  "R2SDHF_21" = "MAD_24_3",
  "R2SDHF_22" = "TX_48_1",
  "R2SDHF_23" = "TX_48_2",
  "R2SDHF_24" = "TX_48_3",
  "R2SDHF_25" = "ROV_48_1",
  "R2SDHF_26" = "ROV_48_2",
  "R2SDHF_27" = "ROV_48_3",
  "R2SDHF_28" = "MAD_48_1",
  "R2SDHF_29" = "MAD_48_2",
  "R2SDHF_30" = "MAD_48_3"
)

# Rename only columns that exist
old <- colnames(counts)
hits <- intersect(old, names(rename_map))
colnames(counts)[match(hits, old)] <- unname(rename_map[hits])

# Checks
setdiff(names(rename_map), colnames(counts))   # should be character(0) if all expected cols present
any(duplicated(colnames(counts)))              # should be FALSE
colnames(counts)


# ---- QC For samples ----
qc_list <- lapply(colnames(counts), function(s) {
  x <- counts[, s, drop = TRUE]
  
  list(
    qc = data.frame(
      sample           = s,
      total_reads      = sum(x),
      detected_genes   = sum(x > 0),
      percent_of_total = sum(x) / sum(colSums(counts)) * 100
    ),
    gene_summary = summary(x)
  )
})
names(qc_list) <- colnames(counts)

qc_summary <- do.call(rbind, lapply(qc_list, `[[`, "qc"))
rownames(qc_summary) <- NULL

gene_summaries <- do.call(cbind, lapply(qc_list, function(z) z$gene_summary))
colnames(gene_summaries) <- names(qc_list)

list(
  qc_table       = qc_summary,
  gene_summaries = gene_summaries
)

# remove low count samples

counts <- counts[, colnames(counts) != "MAD_24_3"]

# ---- set up deseq2 object ----
#set up experimental design
sample_names <- colnames(counts)

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
                       "Mock_0",
                       paste0(group, "_", time))
  )

rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, "Mock_0")



# DESeq model
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = colData,
  design    = ~ condition
)

dds <- DESeq(dds)

resultsNames(dds)

# ---- Shrink data ----
#Look at MA plot
res <- results(dds)
plotMA(res)


#Shrink data

coef_names <- resultsNames(dds)
coef_names <- coef_names[grep("vs_Mock", coef_names)]

res_shrunk_list <- setNames(lapply(coef_names, function(coef_name) {
  lfcShrink(dds, coef = coef_name, type = "apeglm")
}), coef_names)

plotMA(res_shrunk_list[["condition_MAD_12_vs_Mock_0"]])

saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo_resShrink_noMAD24_3.rds")

# ---- PCA ----
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
  Mock_0 = base_cols[["Mock"]],
  TX_12  = tx_shades[1],  TX_24  = tx_shades[2],  TX_48  = tx_shades[3],
  MAD_12 = mad_shades[1], MAD_24 = mad_shades[2], MAD_48 = mad_shades[3],
  ROV_12 = rov_shades[1], ROV_24 = rov_shades[2], ROV_48 = rov_shades[3]
)

pca_data$treatment_time <- factor(
  pca_data$treatment_time,
  levels = c("Mock_0","TX_12","TX_24","TX_48",
             "MAD_12","MAD_24","MAD_48",
             "ROV_12","ROV_24","ROV_48")
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
p

ggsave(
  filename = "figures/plasmo/PCA_treatment_time.png",
  plot = p,
  width = 8,
  height = 6
)

# ---- Volcano Plots ----
padj_cutoff <- 0.05
lfc_cutoff  <- 1.5

outdir <- "figures/plasmo/volcano"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (coef_name in names(res_shrunk_list)) {
  
  volcano_df <- as.data.frame(res_shrunk_list[[coef_name]]) |>
    tibble::rownames_to_column("ensembl_id") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      neglog10_padj = -log10(padj),
      direction = dplyr::case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  n_up   <- sum(volcano_df$direction == "Up", na.rm = TRUE)
  n_down <- sum(volcano_df$direction == "Down", na.rm = TRUE)
  
  x_max <- max(volcano_df$log2FoldChange, na.rm = TRUE)
  x_min <- min(volcano_df$log2FoldChange, na.rm = TRUE)
  y_max <- max(volcano_df$neglog10_padj[is.finite(volcano_df$neglog10_padj)], na.rm = TRUE)
  
  p <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2FoldChange, y = neglog10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = direction), alpha = 0.6, size = 1) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    ggplot2::annotate(
      "text", x = x_max, y = y_max,
      label = paste0("Up: ", n_up),
      hjust = 1, vjust = 1, color = "red"
    ) +
    ggplot2::annotate(
      "text", x = x_min, y = y_max,
      label = paste0("Down: ", n_down),
      hjust = 0, vjust = 1, color = "blue"
    ) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = coef_name,
      x = "log2 Fold Change (apeglm-shrunk)",
      y = "-log10 adjusted p-value"
    )
  
  comp_label <- sub("^condition_", "", coef_name)
  outfile <- file.path(outdir, paste0("volcano_", comp_label, ".png"))
  
  ggplot2::ggsave(outfile, p, width = 7, height = 6, dpi = 300)
}
