#### Joe's Data:Part 2 ####
# Code KG7RYR

# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)
library(readxl)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(tibble)
# ----- Import data ----

#Import expression matrix in tsv format
expr <- read.delim(
  "data/KG7RYR/Plasmo_results/KG7RYR-expression-matrix.tsv",
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
  "KG7RYR_1"  = "control_lung_1",
  "KG7RYR_2"  = "control_lung_2",
  "KG7RYR_3"  = "control_lung_3",
  "KG7RYR_4"  = "control_spleen_1",
  "KG7RYR_5"  = "control_spleen_2",
  "KG7RYR_6"  = "control_spleen_3",
  "KG7RYR_7"  = "cl_lung_4",
  "KG7RYR_8"  = "cl_lung_5",
  "KG7RYR_9"  = "cl_lung_6",
  "KG7RYR_10" = "cl_spleen_4",
  "KG7RYR_11" = "cl_spleen_5",
  "KG7RYR_12" = "cl_spleen_6",
  "KG7RYR_13" = "d_lung_7",
  "KG7RYR_14" = "d_lung_8",
  "KG7RYR_15" = "d_lung_9",
  "KG7RYR_16" = "d_spleen_7",
  "KG7RYR_17" = "d_spleen_8",
  "KG7RYR_18" = "d_spleen_9"
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
counts <- counts[, colnames(counts) != "cl_spleen_6"] # CS12
counts <- counts[, colnames(counts) != "d_spleen_9"] # DS18

#Save/load checkpoint
saveRDS(counts, "data/KG7RYR/r_objects/plasmo_counts(no_CS12_DS18).rds")
counts <- readRDS("data/KG7RYR/r_objects/plasmo_counts(no_CS12_DS18).rds")

# ---- set up deseq2 object ----
# Parse sample metadata from column names
sample_info <- data.frame(
  sample = colnames(counts),
  treatment = gsub("_[a-z]+_[0-9]+$", "", colnames(counts)),
  tissue = gsub("^[a-z]+_([a-z]+)_[0-9]+$", "\\1", colnames(counts)),
  row.names = colnames(counts)
)

# Convert to factors
sample_info$treatment <- factor(sample_info$treatment)
sample_info$tissue <- factor(sample_info$tissue)
sample_info$treatment <- relevel(factor(sample_info$treatment), ref = "control")

# DDs object

dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData = sample_info,
  design = ~ treatment * tissue 
)

# Remove low-information genes to reduce noise and speed up fitting
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# Transform counts for sample-level exploration (PCA)
vsd <- DESeq2::vst(dds)
DESeq2::plotPCA(vsd, intgroup = c("treatment", "tissue")) +
  ggplot2::aes(shape = tissue)

# Fit the DESeq2 model; disable outlier replacement for consistent behavior across conditions
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Shrink all coefficients except the intercept
coef_names <- DESeq2::resultsNames(dds)

coef_names <- DESeq2::resultsNames(dds)[-1]

res_shrunk_list <- setNames(
  lapply(coef_names, function(coef_name) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  }),
  coef_names
)

# Quick check of unshrunk MA plot for the default results
res <- DESeq2::results(dds)
DESeq2::plotMA(res, ylim = c(-5, 5))


#Save/load DESeq object
saveRDS(dds, "data/KG7RYR/r_objects/plasmo_dds_raw.rds")
dds <- readRDS("data/KG7RYR/r_objects/plasmo_dds_raw.rds")

# Normalied reads
dds_norm <- estimateSizeFactors(dds)
normalized_counts <- counts(dds_norm, normalized = TRUE)

saveRDS(normalized_counts, "data/KG7RYR/r_objects/plasmo_counts_normalized.rds")

#Save/load Shrunk Data
saveRDS(res_shrunk_list, "data/KG7RYR/r_objects/plasmo_resShrink_qcmin10.rds")
res_shrunk_list <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_qcmin10.rds")


# ---- ADD mapping of gene names ----
# Ensembl IDs across contrasts ----
all_ensembl <- res_shrunk_list |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Remove version suffix if present (e.g. ENSG00000123456.7)

# ---- Build mapping table ----
gene_map <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = all_ensembl,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)

gene_map <- gene_map |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    symbol     = SYMBOL
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

saveRDS(res_shrunk_list, "data/KG7RYR/r_objects/plasmo_resShrink_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("data/KG7RYR/Plasmo/plasmo_resShrink_qcmin10_annotated.rds")














