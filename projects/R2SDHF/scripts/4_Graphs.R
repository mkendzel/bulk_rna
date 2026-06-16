# ---- Library ----
library(fgsea)
library(ggvenn)
# ---- Data load ----
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)
#

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
  filename = "projects/R2SDHF/figures/plasmo/PCA(nomad)_treatment_time.png",
  plot = p,
  width = 8,
  height = 6
)

# ---- Volcano Plots ----
padj_cutoff <- 0.05
lfc_cutoff  <- 1.5

outdir <- "projects/R2SDHF/figures/plasmo/volcano/nomad"
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

# ---- Heatmap specific pathways ----
library(grid)


res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
# Thresholds used to define “responsive” genes within the selected hallmark set
lfc_cutoff  <- 1.5
padj_cutoff <- 0.05

# Hallmark geneset to visualize
hallmark_name <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
symbols <- unique(hallmark_sets[[hallmark_name]])

# Keep only contrasts that compare to Mock_0
keep_names <- names(res_shrunk_list)[grepl("_vs_Mock_0$", names(res_shrunk_list))]
res_use <- res_shrunk_list[keep_names]

# Normalize Ensembl IDs by stripping version suffix (e.g., ENSG... .12 -> ENSG...)
clean_ens <- function(x) sub("\\..*$", "", as.character(x))

# Detect which columns contain symbols and Ensembl IDs in the DE result tables
df0 <- as.data.frame(res_use[[1]])
sym_col <- intersect(c("SYMBOL", "symbol", "hgnc_symbol"), names(df0))[1]
ens_col <- intersect(c("ensembl_id", "ensembl_gene_id"), names(df0))[1]

# Convert each DE result to a data.frame and keep the contrast name as a column
res_use_df <- lapply(names(res_use), function(nm) {
  df <- as.data.frame(res_use[[nm]])
  df$contrast <- nm
  df
})
names(res_use_df) <- names(res_use)

# Collect Ensembl IDs that are in the hallmark set and significant in any contrast
sig_ens <- unique(unlist(lapply(res_use_df, function(df) {
  df$symbol <- df[[sym_col]]
  df$ensembl_id_clean <- clean_ens(df[[ens_col]])
  df <- df[!is.na(df$symbol) & df$symbol %in% symbols, ]
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df <- df[df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff, ]
  unique(df$ensembl_id_clean)
}), use.names = FALSE))

# Build a stable Ensembl -> symbol mapping from the first contrast table
sym_map <- res_use_df[[1]]
sym_map <- data.frame(
  ensembl_id_clean = clean_ens(sym_map[[ens_col]]),
  symbol = as.character(sym_map[[sym_col]]),
  stringsAsFactors = FALSE
)
sym_map <- sym_map[!is.na(sym_map$ensembl_id_clean) & !duplicated(sym_map$ensembl_id_clean), ]

# Pull DESeq2 normalized counts and add cleaned Ensembl IDs to match DE tables
norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
norm_counts$ensembl_id_clean <- clean_ens(rownames(norm_counts))

# Restrict to genes that are significant in at least one contrast
norm_counts <- norm_counts[norm_counts$ensembl_id_clean %in% sig_ens, , drop = FALSE]
count_cols <- setdiff(names(norm_counts), "ensembl_id_clean")

# Convert wide counts (genes x samples) into long format (one row per gene-sample)
counts_long <- data.frame(
  ensembl_id_clean = rep(norm_counts$ensembl_id_clean, times = length(count_cols)),
  sample = rep(count_cols, each = nrow(norm_counts)),
  norm_count = as.vector(as.matrix(norm_counts[, count_cols, drop = FALSE])),
  stringsAsFactors = FALSE
)

# Bring in sample metadata from the DESeq2 object
coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
coldata_df$sample <- rownames(coldata_df)

# Combine counts, metadata, and gene symbols into a plotting table
plot_df <- dplyr::left_join(counts_long, coldata_df, by = "sample")
plot_df <- dplyr::left_join(plot_df, sym_map, by = "ensembl_id_clean")
plot_df <- plot_df[!is.na(plot_df$symbol) & plot_df$symbol %in% symbols, ]

# Parse the condition label into treatment and time, keeping Mock as time 0
plot_df <- plot_df |>
  dplyr::mutate(
    treatment = dplyr::if_else(grepl("^Mock", condition), "Mock", sub("_.*$", "", condition)),
    time = dplyr::if_else(grepl("^Mock", condition), 0L, as.integer(sub("^.*_", "", condition))),
    group = paste0(treatment, "_", time)
  )

# Set x-axis ordering by time, then by treatment
group_order <- plot_df |>
  dplyr::distinct(group, treatment, time) |>
  dplyr::mutate(
    time = factor(time, levels = sort(unique(time))),
    treatment = factor(treatment, levels = c("Mock", "TX", "ROV"))
  ) |>
  dplyr::arrange(time, treatment) |>
  dplyr::pull(group)

plot_df <- plot_df |>
  dplyr::mutate(group = factor(group, levels = group_order))

# Summarize per gene and group for bar height and error bars
summary_df <- plot_df |>
  dplyr::group_by(symbol, group, treatment) |>
  dplyr::summarise(
    n = sum(!is.na(norm_count)),
    mean_norm = mean(norm_count, na.rm = TRUE),
    sd_norm = stats::sd(norm_count, na.rm = TRUE),
    sem_norm = sd_norm / sqrt(n),
    .groups = "drop"
  )

# Choose whether colors are driven by treatment only or by each group/time combination
color_mode <- "treatment"
file_ext   <- "png"

plot_df <- plot_df |>
  dplyr::mutate(color_key = if (color_mode == "group") as.character(group) else as.character(treatment))

summary_df <- summary_df |>
  dplyr::mutate(color_key = if (color_mode == "group") as.character(group) else as.character(treatment))

# Define fill colors for bars and points
fill_map <- if (color_mode == "treatment") {
  c(Mock = "grey50", TX = "dodgerblue3", ROV = "purple3")
} else {
  setNames(
    rep_len(c("grey50", "dodgerblue3", "purple3", "seagreen3", "orange3"), length(group_order)),
    group_order
  )
}

# Build a bar + jitter plot for a single gene symbol
make_gene_plot <- function(gene_symbol) {
  df_pts <- plot_df |> dplyr::filter(symbol == gene_symbol)
  df_sum <- summary_df |> dplyr::filter(symbol == gene_symbol)
  
  ggplot2::ggplot(df_sum, ggplot2::aes(x = group, y = mean_norm, fill = color_key)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_norm - sem_norm, ymax = mean_norm + sem_norm),
      width = 0.25
    ) +
    ggplot2::geom_point(
      data = df_pts,
      mapping = ggplot2::aes(x = group, y = norm_count, fill = color_key),
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      shape = 21,
      size = 2.5,
      color = "black",
      stroke = 0.6
    ) +
    ggplot2::scale_fill_manual(values = fill_map) +
    ggplot2::labs(x = NULL, y = "Normalized counts", title = gene_symbol) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Gene list to render plots for (one figure per gene)
sig_symbols <- summary_df |>
  dplyr::filter(!is.na(symbol) & symbol != "") |>
  dplyr::distinct(symbol) |>
  dplyr::pull(symbol)

# Output folder per hallmark set
outdir <- file.path("projects/R2SDHF/figures", "bargraphs", hallmark_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Save one plot per gene symbol
for (g in sig_symbols) {
  p <- make_gene_plot(g)
  fn <- file.path(outdir, paste0(g, ".", file_ext))
  ggplot2::ggsave(fn, p, width = 8, height = 4, dpi = 300)
}




# ---- Venn Diagram ----
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
library(grid)
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


library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)

out_dir <- "projects/R2SDHF/data"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()

for (tp in times) {
  
  tx_con  <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TX"]
  rov_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "ROV"]
  
  tx_genes  <- sig_sets[[tx_con]]
  rov_genes <- sig_sets[[rov_con]]
  
  all_genes <- unique(c(tx_genes, rov_genes))
  
  tx_df <- as.data.frame(res_shrunk_list[[tx_con]])
  tx_df$entrez_id <- sub("\\..*$", "", rownames(tx_df))
  tx_df <- tx_df[!duplicated(tx_df$entrez_id), 
                 c("entrez_id", "log2FoldChange", "padj")]
  names(tx_df)[2:3] <- c("TX_log2FC", "TX_padj")
  
  rov_df <- as.data.frame(res_shrunk_list[[rov_con]])
  rov_df$entrez_id <- sub("\\..*$", "", rownames(rov_df))
  rov_df <- rov_df[!duplicated(rov_df$entrez_id), 
                   c("entrez_id", "log2FoldChange", "padj")]
  names(rov_df)[2:3] <- c("ROV_log2FC", "ROV_padj")
  
  out_df <- data.frame(entrez_id = all_genes, stringsAsFactors = FALSE)
  
  out_df$symbol_id <- mapIds(
    org.Hs.eg.db,
    keys      = as.character(out_df$entrez_id),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  
  out_df$membership <- ifelse(
    out_df$entrez_id %in% tx_genes & out_df$entrez_id %in% rov_genes, "Both",
    ifelse(out_df$entrez_id %in% tx_genes, "TX_only", "ROV_only")
  )
  
  out_df <- merge(out_df, tx_df,  by = "entrez_id", all.x = TRUE)
  out_df <- merge(out_df, rov_df, by = "entrez_id", all.x = TRUE)
  
  out_df <- out_df[, c("entrez_id", "symbol_id", "membership",
                       "TX_log2FC", "TX_padj", "ROV_log2FC", "ROV_padj")]
  
  out_df <- out_df[order(factor(out_df$membership, 
                                levels = c("Both", "TX_only", "ROV_only")),
                         out_df$entrez_id), ]
  
  sheet_name <- paste0(tp, "h")
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, out_df)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(out_df)), widths = "auto")
}

saveWorkbook(wb, file.path(out_dir, "DE_genes_venn.xlsx"), overwrite = TRUE)