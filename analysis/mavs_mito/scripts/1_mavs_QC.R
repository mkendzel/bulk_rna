# =============================================
# 1 - QC & data reconciliation
# =============================================
# Reconciles the two raw expression matrices, recovers the missing
# Spleen_D30 samples, runs per-sample QC, drops low-count / old-mock
# samples, and saves the filtered count matrix as `counts_filtered`.
# Output checkpoint feeds 2_mavs_LimmaVoom_DE.R.

# ---- Libraries ----
library(stringr)
library(tidyverse)

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
  lines = c(46:170),
  notes = "Updated raw data to actually include all samples. Specifically, Spleen_D30 samples from 527_EU_Suthar_60RNAseq..."
)


# ---- set up data frame ----
# Identify count columns

expr <- load_checkpoint(
  "all_samples_updated",
  dir = "projects/mavs_mito/data/r_objects"
)


rename_map <- c(
  "6_X_Spleen_D30" = "6_Spleen_D30",
  "18_Y_Spleen_D30" = "18_Spleen_D30",
  "30_Z_Spleen_D30" = "30_Spleen_D30",
  "42_A_Spleen_D30" = "42_Spleen_D30",
  "54_B_Spleen_D30" = "54_Spleen_D30"
  # etc.
)

# Rename only columns that exist
old <- colnames(expr)
hits <- intersect(old, names(rename_map))
colnames(expr)[match(hits, old)] <- unname(rename_map[hits])

# make counts object minus first column of expr
counts <- expr[, -1]
counts[is.na(counts)] <- 0 # replace NA values with 0
unique(is.na(counts$`1_WT_d0`)) # check for NA values in the renamed column

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
counts <- counts[, colnames(counts) != "6_Spleen_D30"]
counts <- counts[, colnames(counts) != "42_Spleen_D30"]
counts <- counts[, colnames(counts) != "18_Spleen_D30"]

#drop old mock
counts <- counts[, colnames(counts) != "1_Spleen_D0"]
counts <- counts[, colnames(counts) != "13_Spleen_D0"]
counts <- counts[, colnames(counts) != "25_Spleen_D0"]
counts <- counts[, colnames(counts) != "37_Spleen_D0"]
counts <- counts[, colnames(counts) != "49_Spleen_D0"]

counts <- counts[, colnames(counts) != "10_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "22_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "34_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "46_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "58_SpleenMAVSko_D0"]

# Save/load checkpoint
save_checkpoint(counts, "counts_filtered", dir = "projects/mavs_mito/data/r_objects",
                lines = c(281:322),
                notes = "Filtered out Spleen samples 6, 18, and 42 at day 30 due to low counts and removed old mock"
)
