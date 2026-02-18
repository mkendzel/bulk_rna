# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)

# ----- Import data ----
#Import expression matrix in tsv format
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
counts <- counts[, colnames(counts) != "MAD_12_1"]
counts <- counts[, colnames(counts) != "MAD_12_2"]
counts <- counts[, colnames(counts) != "MAD_12_3"]
counts <- counts[, colnames(counts) != "MAD_48_1"]
counts <- counts[, colnames(counts) != "MAD_48_2"]
counts <- counts[, colnames(counts) != "MAD_48_3"]
counts <- counts[, colnames(counts) != "MAD_24_1"]
counts <- counts[, colnames(counts) != "MAD_24_2"]
counts <- counts[, colnames(counts) != "MAD_24_3"]

#Save/load checkpoint
# saveRDS(counts, "data/R2SDHF/r_objects/plasmo_counts_noMAD.rds")
# counts <- readRDS("data/R2SDHF/r_objects/plasmo_counts_noMAD.rds")

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

#Save DESeq object
# saveRDS(dds, "data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")
# dds <- readRDS("data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")

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

saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo_resShrink_noMAD.rds")

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
  filename = "figures/plasmo/PCA(nomad)_treatment_time.png",
  plot = p,
  width = 8,
  height = 6
)

# ---- Volcano Plots ----
padj_cutoff <- 0.05
lfc_cutoff  <- 1.5

outdir <- "figures/plasmo/volcano/nomad"
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

# ---- Heatmaps highly variable ----

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

# Ensembl -> SYMBOL map via org.Hs.eg.db
ens_ids <- rownames(mat)

map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL")
) |>
  tibble::as_tibble() |>
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(ensembl_id = ENSEMBL, gene_name = SYMBOL)

idx <- match(ens_ids, map_df$ensembl_id)
mapped_names <- map_df$gene_name[idx]

# Top 50 most variable genes
gene_vars <- apply(mat, 1, var, na.rm = TRUE)
top_ens <- names(sort(gene_vars, decreasing = TRUE))[1:50]

mat_subset <- mat[top_ens, , drop = FALSE]

# Replace Ensembl rownames with gene_name (fallback to Ensembl)
new_names <- mapped_names[match(top_ens, ens_ids)]
new_names[is.na(new_names) | new_names == ""] <- top_ens[is.na(new_names) | new_names == ""]
rownames(mat_subset) <- new_names

# Z-score per gene
mat_scaled <- t(scale(t(mat_subset)))

# Column annotations from colData(dds)

cd <- as.data.frame(colData(dds))

if (!"sample" %in% colnames(cd)) {
  cd$sample <- rownames(cd)
}

annotation_col <- cd |>
  dplyr::mutate(
    condition = as.character(.data$condition),
    treatment = dplyr::if_else(stringr::str_detect(.data$condition, "^Mock"), "Mock",
                               stringr::str_extract(.data$condition, "^[A-Za-z]+")),
    time = dplyr::if_else(.data$treatment == "Mock", "0",
                          stringr::str_extract(.data$condition, "(?<=_)[0-9]+"))
  ) |>
  dplyr::select(sample, treatment, time)

rownames(annotation_col) <- annotation_col$sample
annotation_col$sample <- NULL


# Order columns: Mock first, then TX, MAD, ROV by time (only levels present)
treat_levels <- intersect(c("Mock","TX","ROV"), unique(annotation_col$treatment))

desired_order <- annotation_col |>
  dplyr::mutate(
    treatment = factor(.data$treatment, levels = treat_levels),
    time = as.numeric(.data$time)
  ) |>
  dplyr::arrange(.data$treatment, .data$time) |>
  rownames()

mat_scaled <- mat_scaled[, desired_order, drop = FALSE]
annotation_col <- annotation_col[desired_order, , drop = FALSE]

pheatmap::pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE
)



# ---- log2FC heatmap vs Mock_0 (non-sig black), biomaRt gene labels, no MAD ----
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

lfc_cutoff  <- 1.5
padj_cutoff <- 0.05
top_n_genes <- 250

# Identify contrasts comparing each condition to Mock_0
keep_names <- names(res_shrunk_list)[stringr::str_detect(names(res_shrunk_list), "_vs_Mock_0$")]
res_use <- res_shrunk_list[keep_names]

# Convert DESeqResults objects to data frames and retain contrast names
res_use_df <- purrr::imap(res_use, function(res, contrast_name) {
  df <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  df$contrast   <- contrast_name
  df
})

# Combine all contrasts and retain genes meeting significance thresholds
sig_tbl <- dplyr::bind_rows(res_use_df) |>
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
  dplyr::filter(padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) |>
  dplyr::select(ensembl_id, log2FoldChange, padj, contrast)

# Rank genes by maximum absolute log2 fold change across contrasts
top_genes <- sig_tbl |>
  dplyr::group_by(ensembl_id) |>
  dplyr::summarise(max_abs_lfc = max(abs(log2FoldChange)), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(max_abs_lfc)) |>
  dplyr::slice_head(n = top_n_genes) |>
  dplyr::pull(ensembl_id)

# Construct matrix of log2 fold changes for selected genes
# Values not meeting thresholds are set to NA for visualization
lfc_mat <- sapply(res_use_df, function(df) {
  idx <- match(top_genes, df$ensembl_id)
  lfc <- df$log2FoldChange[idx]
  p   <- df$padj[idx]
  
  is_sig <- !is.na(p) & !is.na(lfc) & (p < padj_cutoff) & (abs(lfc) >= lfc_cutoff)
  lfc[!is_sig] <- NA_real_
  
  lfc
})
colnames(lfc_mat) <- keep_names

# Parse contrast names and generate simplified column labels (e.g., "TX 12")
contrast_tbl <- tibble::tibble(contrast = colnames(lfc_mat)) |>
  dplyr::mutate(
    condition = stringr::str_match(contrast, "condition_([^_]+_[0-9]+)_vs_Mock_0")[, 2],
    treatment = stringr::str_extract(condition, "^[A-Za-z]+"),
    time      = as.numeric(stringr::str_extract(condition, "(?<=_)[0-9]+")),
    label     = paste(treatment, time)
  )

# Order columns by time, then treatment (TX and ROV adjacent within each timepoint)
treat_levels <- c("TX", "ROV")
contrast_order <- contrast_tbl |>
  dplyr::mutate(
    time = factor(time, levels = sort(unique(time))),
    treatment = factor(treatment, levels = c(treat_levels, setdiff(unique(treatment), treat_levels)))
  ) |>
  dplyr::arrange(time, treatment) |>
  dplyr::pull(contrast)

lfc_mat <- lfc_mat[, contrast_order, drop = FALSE]

# Apply simplified column labels
col_label_map <- contrast_tbl$label
names(col_label_map) <- contrast_tbl$contrast
colnames(lfc_mat) <- col_label_map[colnames(lfc_mat)]

# Map Ensembl identifiers to HGNC symbols; fallback to original ID if unavailable
top_genes_clean <- sub("\\..*$", "", top_genes)

bm <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = unique(top_genes_clean),
  mart       = mart
) |>
  tibble::as_tibble() |>
  dplyr::mutate(gene_name = dplyr::coalesce(hgnc_symbol, external_gene_name)) |>
  dplyr::filter(!is.na(gene_name) & gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

gene_labels <- bm$gene_name[match(top_genes_clean, bm$ensembl_gene_id)]
gene_labels[is.na(gene_labels) | gene_labels == ""] <- top_genes[is.na(gene_labels) | gene_labels == ""]
rownames(lfc_mat) <- gene_labels

# Replace NA with zero only for clustering distance calculation
lfc_for_cluster <- lfc_mat
lfc_for_cluster[is.na(lfc_for_cluster)] <- 0
row_hc <- stats::hclust(stats::dist(lfc_for_cluster))

# Generate heatmap with non-significant values displayed in black
pheatmap::pheatmap(
  lfc_mat,
  cluster_rows = row_hc,
  cluster_cols = FALSE,
  na_col = "black",
  show_colnames = TRUE,
  show_rownames = FALSE,
  angle_col = 0,
  main = "Top 250 DE Genes Relative to Mock (log2 Fold Change)"
)

# ---- Venn Diagram of DE Genes ----

lfc_cutoff  <- 1.5
padj_cutoff <- 0.05

# Identify contrasts vs Mock_0
keep_names <- names(res_shrunk_list)
keep_names <- keep_names[grepl("_vs_Mock_0$", keep_names)]

# Parse treatment and time from contrast names
contrast_tbl <- data.frame(
  contrast  = keep_names,
  condition = sub("condition_([^_]+_[0-9]+)_vs_Mock_0", "\\1", keep_names),
  stringsAsFactors = FALSE
)

contrast_tbl$treatment <- sub("_.*", "", contrast_tbl$condition)
contrast_tbl$time      <- as.numeric(sub(".*_", "", contrast_tbl$condition))
contrast_tbl <- contrast_tbl[contrast_tbl$treatment %in% c("TX", "ROV"), ]

# Build significant gene sets per contrast
sig_sets <- list()

for (nm in keep_names) {
  res <- res_shrunk_list[[nm]]
  df  <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  
  sig_genes <- df$ensembl_id[
    !is.na(df$padj) &
      !is.na(df$log2FoldChange) &
      df$padj < padj_cutoff &
      abs(df$log2FoldChange) >= lfc_cutoff
  ]
  
  sig_sets[[nm]] <- unique(sub("\\..*$", "", sig_genes))
}

# Create Venn plots for 12, 24, 48 hours
times <- sort(unique(contrast_tbl$time))
times <- times[times %in% c(12, 24, 48)]

venn_plots <- list()

for (tp in times) {
  
  tx_con  <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TX"]
  rov_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "ROV"]
  
  sets_tp <- list(
    TX  = sig_sets[[tx_con]],
    ROV = sig_sets[[rov_con]]
  )
  
  venn_plots[[as.character(tp)]] <-
    ggvenn(
      sets_tp,
      show_elements = FALSE
    ) +
    ggtitle(paste0(tp, "h")) +
    theme(plot.title = element_text(hjust = 0.5),
          size = 40)
}

# Draw all three in a single row
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(venn_plots))))

for (i in seq_along(venn_plots)) {
  print(
    venn_plots[[i]],
    vp = viewport(layout.pos.row = 1, layout.pos.col = i)
  )
}
