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

# =====================================
# ---- set up deseq2 object (FULL) ----
# =====================================
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


# ---- Annotate genes to symbol id ----
all_ensembl <- res_shrunk_list |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Build mapping table
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

# Append identifiers to every contrast
res_shrunk_list <- purrr::map(res_shrunk_list, function(df) {
  
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map, by = "ensembl_id") |>
    as_tibble()
  
})

# Save and load annotated list
saveRDS(res_shrunk_list, "data/KG7RYR/r_objects/plasmo_resShrink_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("data/KG7RYR/Plasmo/plasmo_resShrink_qcmin10_annotated.rds")


# ==================================
# ---- Tissue-specific analyses ----
# ==================================
# ---- Lung ----

#Subset to lung
lung_samples <- sample_info$tissue == "lung"
counts_lung <- counts[, lung_samples]
info_lung <- sample_info[lung_samples, ]
info_lung$treatment <- relevel(droplevels(info_lung$treatment), ref = "control")

# DESEQ object creation
dds_lung <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_lung)),
  colData = info_lung,
  design = ~ treatment
)

#Subset low expression genes out
keep <- rowSums(DESeq2::counts(dds_lung) >= 10) >= 5
dds_lung <- dds_lung[keep, ]

#variance stabilization of pca
vsd_lung <- DESeq2::vst(dds_lung)
DESeq2::plotPCA(vsd_lung, intgroup = "treatment")

#run DESeq2
dds_lung <- DESeq2::DESeq(dds_lung, minReplicatesForReplace = Inf)

# Shrink logfold2 estimations
coef_names_lung <- DESeq2::resultsNames(dds_lung)[-1]
res_shrunk_lung <- setNames(
  lapply(coef_names_lung, function(coef_name) {
    DESeq2::lfcShrink(dds_lung, coef = coef_name, type = "apeglm")
  }),
  coef_names_lung
)

#Add all comps to results
res_wald_lung <- setNames(
  lapply(coef_names_lung, function(coef_name) {
    DESeq2::results(dds_lung, name = coef_name)
  }),
  coef_names_lung
)
res_wald_lung[["treatment_cl_vs_d"]] <- DESeq2::results(dds_lung, contrast = c("treatment", "cl", "d"))

res_lung <- DESeq2::results(dds_lung)
DESeq2::plotMA(res_lung, ylim = c(-5, 5))

#Save/load DESeq object
saveRDS(dds_lung, "data/KG7RYR/r_objects/plasmo_dds_lung_raw.rds")
dds_lung <- readRDS("data/KG7RYR/r_objects/plasmo_dds_lung_raw.rds")

saveRDS(counts(dds_lung, normalized = TRUE), "data/KG7RYR/r_objects/plasmo_counts_lung_normalized.rds")
normalized_counts_lung <- readRDS("data/KG7RYR/r_objects/plasmo_counts_lung_normalized.rds")

saveRDS(res_shrunk_lung, "data/KG7RYR/r_objects/plasmo_resShrink_lung_qcmin10.rds")
res_shrunk_lung <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_lung_qcmin10.rds")

saveRDS(res_wald_lung, "data/KG7RYR/r_objects/plasmo_resWald_lung_qcmin10.rds")
res_wald_lung <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_lung_qcmin10.rds")

# Annotated object
all_ensembl_lung <- res_shrunk_lung |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

gene_map_lung <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = all_ensembl_lung,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
gene_map_lung <- gene_map_lung |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    symbol     = SYMBOL
  )

res_shrunk_lung <- purrr::map(res_shrunk_lung, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_lung, by = "ensembl_id") |>
    as_tibble()
})

res_wald_lung <- purrr::map(res_wald_lung, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_lung, by = "ensembl_id") |>
    as_tibble()
})

# Save Annotated
saveRDS(res_shrunk_lung, "data/KG7RYR/r_objects/plasmo_resShrink_lung_qcmin10_annotated.rds")
res_shrunk_lung <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_lung_qcmin10_annotated.rds")

saveRDS(res_wald_lung, "data/KG7RYR/r_objects/plasmo_resWald_lung_qcmin10_annotated.rds")
res_wald_lung <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_lung_qcmin10_annotated.rds")

# ---- Spleen ----
# Subsample to spleen
spleen_samples <- sample_info$tissue == "spleen"
counts_spleen <- counts[, spleen_samples]
info_spleen <- sample_info[spleen_samples, ]
info_spleen$treatment <- relevel(droplevels(info_spleen$treatment), ref = "control")

# Create DESEQ2 object
dds_spleen <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_spleen)),
  colData = info_spleen,
  design = ~ treatment
)

# Remove lowly expressed genes
keep <- rowSums(DESeq2::counts(dds_spleen) >= 10) >= 5
dds_spleen <- dds_spleen[keep, ]

#Variance stabalization for PCA
vsd_spleen <- DESeq2::vst(dds_spleen)
DESeq2::plotPCA(vsd_spleen, intgroup = "treatment")

# Run DESeq2
dds_spleen <- DESeq2::DESeq(dds_spleen, minReplicatesForReplace = Inf)

#Shrunken results
coef_names_spleen <- DESeq2::resultsNames(dds_spleen)[-1]
res_shrunk_spleen <- setNames(
  lapply(coef_names_spleen, function(coef_name) {
    DESeq2::lfcShrink(dds_spleen, coef = coef_name, type = "apeglm")
  }),
  coef_names_spleen
)

#Wald test stat results (full comp)
res_wald_spleen <- setNames(
  lapply(coef_names_spleen, function(coef_name) {
    DESeq2::results(dds_spleen, name = coef_name)
  }),
  coef_names_spleen
)
res_wald_spleen[["treatment_cl_vs_d"]] <- DESeq2::results(dds_spleen, contrast = c("treatment", "cl", "d"))

res_spleen <- DESeq2::results(dds_spleen)
DESeq2::plotMA(res_spleen, ylim = c(-5, 5))

#Save and load checkpoints
saveRDS(dds_spleen, "data/KG7RYR/r_objects/plasmo_dds_spleen_raw.rds")
dds_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_dds_spleen_raw.rds")

saveRDS(counts(dds_spleen, normalized = TRUE), "data/KG7RYR/r_objects/plasmo_counts_spleen_normalized.rds")
normalized_counts_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_counts_spleen_normalized.rds")

saveRDS(res_shrunk_spleen, "data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10.rds")
res_shrunk_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10.rds")

saveRDS(res_wald_spleen, "data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10.rds")
res_wald_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10.rds")

#Annotated

all_ensembl_spleen <- res_shrunk_spleen |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

gene_map_spleen <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = all_ensembl_spleen,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
gene_map_spleen <- gene_map_spleen |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    symbol     = SYMBOL
  )

res_shrunk_spleen <- purrr::map(res_shrunk_spleen, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_spleen, by = "ensembl_id") |>
    as_tibble()
})

res_wald_spleen <- purrr::map(res_wald_spleen, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_spleen, by = "ensembl_id") |>
    as_tibble()
})

# save/load
saveRDS(res_shrunk_spleen, "data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")
res_shrunk_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")

saveRDS(res_wald_spleen, "data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")
res_wald_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")














