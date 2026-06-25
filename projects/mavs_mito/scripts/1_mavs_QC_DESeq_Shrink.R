# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)
library(readxl)
# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

# ----- Import data ----
## Primary import: expression matrix (.txt)
expr <- read.delim(
  "projects/mavs_mito/data/raw/All samples.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(expr) <- expr$gene_id
expr <- expr[, -1]

## Alternative import: subset of expression matrix (.txt)
sub_expr <- read.delim(
  "projects/mavs_mito/data/raw/527_EU_Suthar_60RNAseq_4_RawCounts_wGeneSymbols.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#remove first four rows
sub_expr <- sub_expr[-(1:4), ]

#Rename rownames to gene_id
rownames(sub_expr) <- sub_expr$gene_id
sub_expr <- sub_expr[, -1]


# ---- Explore samples ----
## 1) Count unique samples in each data set
# "all data set"
get_suffix <- function(cols) sub("^[0-9]+_", "", cols)

expr_suffixes <- get_suffix(colnames(expr))

suffix_counts <- as.data.frame(table(expr_suffixes), stringsAsFactors = FALSE)
names(suffix_counts) <- c("suffix", "n")
suffix_counts <- suffix_counts[order(-suffix_counts$n), ]
rownames(suffix_counts) <- NULL

suffix_counts

# "subset data set"
get_suffix_drop2 <- function(cols) sub("^[^_]*_[^_]*_", "", cols)

sub_expr_suffixes <- get_suffix_drop2(colnames(sub_expr))

sub_suffix_counts <- as.data.frame(table(sub_expr_suffixes), stringsAsFactors = FALSE)
names(sub_suffix_counts) <- c("suffix", "n")
sub_suffix_counts <- sub_suffix_counts[order(-sub_suffix_counts$n), ]
rownames(sub_suffix_counts) <- NULL

sub_suffix_counts

# Table of suffixes in both data sets
combined_counts <- merge(
  suffix_counts, sub_suffix_counts,
  by = "suffix", all = TRUE,
  suffixes = c("_expr", "_sub")
)
names(combined_counts) <- c("suffix", "n_expr", "n_sub")
combined_counts[is.na(combined_counts)] <- 0
combined_counts <- combined_counts[order(-combined_counts$n_expr, -combined_counts$n_sub), ]
rownames(combined_counts) <- NULL

knitr::kable(combined_counts, format = "markdown")

# Notes:
# Missing 3 Spleen Day 30 samples in the "All data set" stored in the object "expr"
# Solution: Include the missing samples to "All data" from the "subset data set"
# Pass these through QC to see why they were removed

# ---- Compare genes and gene counts between data sets ----
# Extract leading sample number from column names
get_sample_num <- function(cols) as.integer(sub("^([0-9]+).*", "\\1", cols))

expr_num <- get_sample_num(colnames(expr))
sub_num <- get_sample_num(colnames(sub_expr))

expr_lookup <- data.frame(col = colnames(expr), num = expr_num)
expr_lookup <- expr_lookup[!is.na(expr_lookup$num), ]

sub_lookup <- data.frame(col = colnames(sub_expr), num = sub_num)
sub_lookup <- sub_lookup[!is.na(sub_lookup$num), ]

## 1) Do the two data sets share the same samples? ----
samples_only_in_expr <- sort(setdiff(expr_lookup$num, sub_lookup$num))
samples_only_in_sub <- sort(setdiff(sub_lookup$num, expr_lookup$num))
samples_dup_in_expr <- sort(unique(expr_lookup$num[duplicated(expr_lookup$num)]))
samples_dup_in_sub <- sort(unique(sub_lookup$num[duplicated(sub_lookup$num)]))

same_samples <- length(samples_only_in_expr) == 0 &&
  length(samples_only_in_sub) == 0 &&
  length(samples_dup_in_expr) == 0 &&
  length(samples_dup_in_sub) == 0

sample_check <- list(
  same_samples = same_samples,
  samples_only_in_expr = samples_only_in_expr,
  samples_only_in_sub = samples_only_in_sub,
  duplicated_sample_numbers_in_expr = samples_dup_in_expr,
  duplicated_sample_numbers_in_sub = samples_dup_in_sub
)
sample_check

# Sample numbers safe to compare: present in both, not duplicated in either
shared_nums <- sort(setdiff(
  intersect(expr_lookup$num, sub_lookup$num),
  union(samples_dup_in_expr, samples_dup_in_sub)
))

## 2) Within shared samples, how many genes are shared? ----
shared_genes <- intersect(rownames(expr), rownames(sub_expr))

gene_check <- list(
  n_genes_expr = nrow(expr),
  n_genes_sub = nrow(sub_expr),
  n_genes_shared = length(shared_genes),
  genes_only_in_expr = setdiff(rownames(expr), rownames(sub_expr)),
  genes_only_in_sub = setdiff(rownames(sub_expr), rownames(expr))
)
gene_check

## 3) Do shared genes have identical counts? ----
# not present as comparable conditions in sub_expr
cols_to_remove <- c(
  "1_WT_d0", "2_WT_d0", "3_WT_d0", "4_WT_d0",
  "5_MAVSKO_d0", "6_MAVSKO_d0", "7_MAVSKO_d0", "8_MAVSKO_d0"
)

expr_filtered <- expr[, !(colnames(expr) %in% cols_to_remove)]

# Rebuild sample lookup tables using the filtered expr
expr_num <- get_sample_num(colnames(expr_filtered))
expr_lookup <- data.frame(col = colnames(expr_filtered), num = expr_num)
expr_lookup <- expr_lookup[!is.na(expr_lookup$num), ]

sub_num <- get_sample_num(colnames(sub_expr))
sub_lookup <- data.frame(col = colnames(sub_expr), num = sub_num)
sub_lookup <- sub_lookup[!is.na(sub_lookup$num), ]

# Recheck for duplicates/unmatched samples after filtering
samples_only_in_expr <- sort(setdiff(expr_lookup$num, sub_lookup$num))
samples_only_in_sub <- sort(setdiff(sub_lookup$num, expr_lookup$num))
samples_dup_in_expr <- sort(unique(expr_lookup$num[duplicated(expr_lookup$num)]))
samples_dup_in_sub <- sort(unique(sub_lookup$num[duplicated(sub_lookup$num)]))

shared_nums <- sort(setdiff(
  intersect(expr_lookup$num, sub_lookup$num),
  union(samples_dup_in_expr, samples_dup_in_sub)
))

# Shared genes between filtered expr and sub_expr
shared_genes <- intersect(rownames(expr_filtered), rownames(sub_expr))

compare_one_sample <- function(n) {
  expr_col <- expr_lookup$col[expr_lookup$num == n]
  sub_col <- sub_lookup$col[sub_lookup$num == n]

  x <- as.numeric(expr_filtered[shared_genes, expr_col])
  y <- as.numeric(sub_expr[shared_genes, sub_col])
  is_match <- x == y

  data.frame(
    number = n,
    expr_col = expr_col,
    sub_col = sub_col,
    n_genes = length(shared_genes),
    n_mismatch = sum(!is_match, na.rm = TRUE),
    pct_mismatch = round(mean(!is_match, na.rm = TRUE) * 100, 2)
  )
}

sample_comparison <- do.call(rbind, lapply(shared_nums, compare_one_sample))
sample_comparison <- sample_comparison[order(sample_comparison$number), ]
rownames(sample_comparison) <- NULL

sample_comparison
mismatched_samples <- sample_comparison[sample_comparison$n_mismatch > 0, ]
mismatched_samples

# Notes: There are no missmatched counts. The data sets are identical except for the missing Spleen day 30s

# ---- Edit expr to include missing samples ----
# Remove all Spleen_D30 columns from expr
expr_no_d30 <- expr[, !grepl("Spleen_D30", colnames(expr))]

# Add all Spleen_D30 columns from sub_expr
d30_cols_sub <- colnames(sub_expr)[grepl("Spleen_D30", colnames(sub_expr))]
expr_updated <- cbind(expr_no_d30, sub_expr[rownames(expr_no_d30), d30_cols_sub])

#Check
get_suffix <- function(cols) sub("^[0-9]+_", "", cols)

expr_suffixes <- get_suffix(colnames(expr_updated))

suffix_counts <- as.data.frame(table(expr_suffixes), stringsAsFactors = FALSE)
names(suffix_counts) <- c("suffix", "n")
suffix_counts <- suffix_counts[order(-suffix_counts$n), ]
rownames(suffix_counts) <- NULL

suffix_counts 

# Save checkpoint
save_checkpoint(
  expr_updated,
  "all_samples_updated",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(12:99),
  notes = "Updated raw data to actually include all samples. Specifically, Spleen_D30 samples from 527_EU_Suthar_60RNAseq..."
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


# rename columns to project-specific names.

rename_map <- c(
  "1_Spleen_D0" = "spleen_mock",
  "2_Spleen_D7" = "spleen_d7",
  "3_Brain_D7" = "brain_d7",
  "4_Spleen_D15" = "spleen_d15",
  "5_Brain_D15" = "brain_d15",
  ""
  # etc.
)

# Rename only columns that exist
old <- colnames(counts)
hits <- intersect(old, names(rename_map))
colnames(counts)[match(hits, old)] <- unname(rename_map[hits])

# Checks
setdiff(names(rename_map), colnames(counts)) # should be character(0) if all expected cols present
any(duplicated(colnames(counts))) # should be FALSE
colnames(counts)


# ---- QC For samples ----
qc_list <- lapply(colnames(counts), function(s) {
  x <- counts[, s, drop = TRUE]

  list(
    qc = data.frame(
      sample = s,
      total_reads = sum(x),
      detected_genes = sum(x > 0),
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
  qc_table = qc_summary,
  gene_summaries = gene_summaries
)

# Drop low-count samples based on qc_summary/qc_list above
counts <- counts[, colnames(counts) != "TreatmentB_1_1"]
counts <- counts[, colnames(counts) != "TreatmentB_1_2"]


# Save/load checkpoint
# save_checkpoint(counts, "counts_filtered", dir = "projects/YOUR_PROJECT_NAME/data/r_objects")
counts <- load_checkpoint(
  "counts_filtered",
  dir = "projects/YOUR_PROJECT_NAME/data/r_objects"
)

# ---- set up deseq2 object ----
sample_names <- colnames(counts)

# Strip trailing replicate number to get condition (e.g. TreatmentA_1_2 -> TreatmentA_1)
colData <- data.frame(sample = sample_names) |>
  dplyr::mutate(
    condition = sub("_[0-9]+$", "", sample)
  )

rownames(colData) <- colData$sample

# EDIT: set your control/reference condition
control_condition <- "Control_0"
colData$condition <- relevel(factor(colData$condition), ref = control_condition)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = colData,
  design = ~condition
)

# Filter low-count genes with edgeR
keep <- edgeR::filterByExpr(dds, group = colData$condition)
dds <- dds[keep, ]

# Transform counts for sample-level exploration (PCA)
vsd <- DESeq2::vst(dds)
DESeq2::plotPCA(vsd, intgroup = "condition")

# Fit the model, disable outlier replacement
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Inspect the coefficient names to confirm which contrasts exist
coef_names <- DESeq2::resultsNames(dds)

# Keep only coefficients vs control
coef_names_vs_ref <- coef_names[grep(
  paste0("vs_", control_condition),
  coef_names
)]

# Shrink LFCs with apeglm for each contrast
res_shrunk_list <- setNames(
  lapply(coef_names_vs_ref, function(coef_name) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  }),
  coef_names_vs_ref
)

# Quick check of unshrunk MA plot for the default results
res <- DESeq2::results(dds)
DESeq2::plotMA(res, ylim = c(-5, 5))


# Save/load DESeq object
# save_checkpoint(dds, "dds_filtered", dir = "projects/PROJECT_NAME/data/r_objects")
dds <- load_checkpoint(
  "dds_filtered",
  dir = "projects/PROJECT_NAME/data/r_objects"
)

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

# Save/load normalized counts
save_checkpoint(
  normalized_counts,
  "normalized_counts",
  dir = "projects/PROJECT_NAME/data/r_objects"
)
normalized_counts <- load_checkpoint(
  "normalized_counts",
  dir = "projects/PROJECT_NAME/data/r_objects"
)

# Save/load Shrunk Data
save_checkpoint(
  res_shrunk_list,
  "res_shrunk_list",
  dir = "projects/PROJECT_NAME/data/r_objects"
)
res_shrunk_list <- load_checkpoint(
  "res_shrunk_list",
  dir = "projects/PROJECT_NAME/data/r_objects"
)


# ---- Map gene IDs ----
all_ensembl <- res_shrunk_list |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Remove version suffix if present (e.g. ENSG00000123456.7)
all_ensembl_clean <- sub("\\..*$", "", all_ensembl)

# ---- Build mapping table ----
# Human: org.Hs.eg.db. Mouse: org.Mm.eg.db
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_ensembl_clean,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  rename(
    ensembl_id = ENSEMBL,
    entrez_id = ENTREZID,
    SYMBOL = SYMBOL
  )

# ---- Append identifiers to every contrast ----
res_shrunk_list <- purrr::map(res_shrunk_list, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map, by = "ensembl_id") |>
    as_tibble()
})

save_checkpoint(
  res_shrunk_list,
  "res_shrunk_list_annotated",
  dir = "projects/PROJECT_NAME/data/r_objects"
)
res_shrunk_list <- load_checkpoint(
  "res_shrunk_list_annotated",
  dir = "projects/PROJECT_NAME/data/r_objects"
)
