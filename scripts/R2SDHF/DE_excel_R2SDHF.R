library(tidyverse)
library(DESeq2)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(openxlsx)

# Load DE results and DESeq2 object
res_shrunk_list <- readRDS("data/R2SDHF/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
dds <- readRDS("data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")

# You'll need the dds object loaded - adjust path as needed
# dds <- readRDS("path/to/your/dds.rds")

# Hallmark gene sets
hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Thresholds
lfc_cutoff  <- 1.5
padj_cutoff <- 0.05

# Hallmark geneset
hallmark_name <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
symbols <- unique(hallmark_sets[[hallmark_name]])

# Keep only contrasts vs Mock_0
keep_names <- names(res_shrunk_list)[grepl("_vs_Mock_0$", names(res_shrunk_list))]
res_use <- res_shrunk_list[keep_names]

# Helper to strip Ensembl version suffix
clean_ens <- function(x) sub("\\..*$", "", as.character(x))

# Detect column names for symbol and ensembl ID
df0 <- as.data.frame(res_use[[1]])
sym_col <- intersect(c("SYMBOL", "symbol", "hgnc_symbol"), names(df0))[1]
ens_col <- intersect(c("ensembl_id", "ensembl_gene_id"), names(df0))[1]

# Convert each DE result to data.frame
res_use_df <- lapply(names(res_use), function(nm) {
  df <- as.data.frame(res_use[[nm]])
  df$contrast <- nm
  df
})
names(res_use_df) <- names(res_use)

# Collect Ensembl IDs significant in any contrast
sig_ens <- unique(unlist(lapply(res_use_df, function(df) {
  df$symbol <- df[[sym_col]]
  df$ensembl_id_clean <- clean_ens(df[[ens_col]])
  df <- df[!is.na(df$symbol) & df$symbol %in% symbols, ]
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df <- df[df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff, ]
  unique(df$ensembl_id_clean)
}), use.names = FALSE))

# Ensembl -> symbol mapping
sym_map <- res_use_df[[1]]
sym_map <- data.frame(
  ensembl_id_clean = clean_ens(sym_map[[ens_col]]),
  symbol = as.character(sym_map[[sym_col]]),
  stringsAsFactors = FALSE
)
sym_map <- sym_map[!is.na(sym_map$ensembl_id_clean) & !duplicated(sym_map$ensembl_id_clean), ]

# Normalized counts
norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
norm_counts$ensembl_id_clean <- clean_ens(rownames(norm_counts))
norm_counts <- norm_counts[norm_counts$ensembl_id_clean %in% sig_ens, , drop = FALSE]
count_cols <- setdiff(names(norm_counts), "ensembl_id_clean")

# Pivot to long format
counts_long <- data.frame(
  ensembl_id_clean = rep(norm_counts$ensembl_id_clean, times = length(count_cols)),
  sample = rep(count_cols, each = nrow(norm_counts)),
  norm_count = as.vector(as.matrix(norm_counts[, count_cols, drop = FALSE])),
  stringsAsFactors = FALSE
)

# Merge metadata and symbols
coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
coldata_df$sample <- rownames(coldata_df)

plot_df <- left_join(counts_long, coldata_df, by = "sample") %>%
  left_join(sym_map, by = "ensembl_id_clean") %>%
  filter(!is.na(symbol), symbol %in% symbols) %>%
  mutate(
    treatment = if_else(grepl("^Mock", condition), "Mock", sub("_.*$", "", condition)),
    time = if_else(grepl("^Mock", condition), 0L, as.integer(sub("^.*_", "", condition))),
    group = paste0(treatment, "_", time)
  )

# Set group ordering
group_order <- plot_df %>%
  distinct(group, treatment, time) %>%
  mutate(
    time = factor(time, levels = sort(unique(time))),
    treatment = factor(treatment, levels = c("Mock", "TX", "ROV"))
  ) %>%
  arrange(time, treatment) %>%
  pull(group)

plot_df <- plot_df %>%
  mutate(group = factor(group, levels = group_order))

# Summary stats
summary_df <- plot_df %>%
  group_by(symbol, group, treatment) %>%
  summarise(
    n = sum(!is.na(norm_count)),
    mean_norm = mean(norm_count, na.rm = TRUE),
    sd_norm = sd(norm_count, na.rm = TRUE),
    sem_norm = sd_norm / sqrt(n),
    .groups = "drop"
  )

# Gene list
sig_symbols <- summary_df %>%
  filter(!is.na(symbol), symbol != "") %>%
  distinct(symbol) %>%
  pull(symbol)

# Build Excel workbook with one tab per gene
wb <- createWorkbook()

for (g in sig_symbols) {
  df_pts <- plot_df %>%
    filter(symbol == g) %>%
    dplyr::select(sample, group, treatment, time, norm_count)
  
  df_sum <- summary_df %>%
    filter(symbol == g) %>%
    dplyr::select(group, treatment, n, mean_norm, sd_norm, sem_norm)
  
  sheet_name <- substr(gsub("[^A-Za-z0-9]", "_", g), 1, 31)
  
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, "Individual Sample Counts", startRow = 1)
  writeData(wb, sheet_name, df_pts, startRow = 2)
  
  gap_row <- nrow(df_pts) + 4
  writeData(wb, sheet_name, "Summary Statistics", startRow = gap_row)
  writeData(wb, sheet_name, df_sum, startRow = gap_row + 1)
}

outdir <- file.path("figures", "bargraphs", hallmark_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
saveWorkbook(wb, file.path(outdir, paste0(hallmark_name, "_counts.xlsx")), overwrite = TRUE)
