# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)
library(readxl)
# ----- Import data ----
## Shilu's data
#Import expression matrix in tsv format
expr <- read.delim(
  "projects/R2SDHF/data/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## Jimina's data
#Import 
de_files <- list.files("projects/Anderson-Suthar_collab/data", pattern = "^DE_", full.names = TRUE)

de_list <- setNames(
  lapply(de_files, function(f) {
    if (grepl("\\.csv$", f)) read.csv(f) else read_excel(f)
  }),
  gsub("\\.(xlsx|csv)$", "", basename(de_files))
)

names(de_list)



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
# saveRDS(counts, "projects/R2SDHF/data/r_objects/plasmo_counts_noMAD.rds")
counts <- readRDS("projects/R2SDHF/data/r_objects/plasmo_counts_noMAD.rds")

# ---- set up deseq2 object ----
sample_names <- colnames(counts)

colData <- data.frame(sample = sample_names) |>
  dplyr::mutate(
    group = dplyr::case_when(
      stringr::str_detect(sample, "Mock") ~ "Mock",
      TRUE ~ stringr::str_extract(sample, "TX|MAD|ROV")
    ),
    time = dplyr::case_when(
      stringr::str_detect(sample, "Mock") ~ "0",
      TRUE ~ stringr::str_extract(sample, "12|24|48")
    ),
    condition = ifelse(group == "Mock", "Mock_0", paste0(group, "_", time))
  )

rownames(colData) <- colData$sample

# Ensure the reference level is Mock_0 so all contrasts are computed vs Mock_0
colData$condition <- relevel(factor(colData$condition), ref = "Mock_0")

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = colData,
  design    = ~ condition
)

# Remove low-information genes to reduce noise and speed up fitting
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# Transform counts for sample-level exploration (PCA)
vsd <- DESeq2::vst(dds)
DESeq2::plotPCA(vsd, intgroup = "condition")

# Fit the DESeq2 model; disable outlier replacement for consistent behavior across conditions
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Inspect the coefficient names to confirm which contrasts exist
coef_names <- DESeq2::resultsNames(dds)

# Keep only coefficients that represent contrasts against the Mock_0 reference
coef_names_vs_mock <- coef_names[grep("vs_Mock", coef_names)]

# Shrink log2 fold-changes with apeglm for each vs-Mock coefficient
res_shrunk_list <- setNames(
  lapply(coef_names_vs_mock, function(coef_name) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  }),
  coef_names_vs_mock
)

# Quick check of unshrunk MA plot for the default results
res <- DESeq2::results(dds)
DESeq2::plotMA(res, ylim = c(-5, 5))


#Save/load DESeq object
# saveRDS(dds, "projects/R2SDHF/data/r_objects/plasmo_dds_noMAD.rds")
dds <- readRDS("projects/R2SDHF/data/r_objects/plasmo_dds_noMAD.rds")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#Save/load Shrunk Data
saveRDS(res_shrunk_list, "projects/R2SDHF/data/Plasmo_resShrink_noMAD_qcmin10.rds")
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo_resShrink_noMAD_qcmin10.rds")


# ---- ADD mapping of gene names ----
# Ensembl IDs across contrasts ----
all_ensembl <- res_shrunk_list |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Remove version suffix if present (e.g. ENSG00000123456.7)
all_ensembl_clean <- sub("\\..*$", "", all_ensembl)

# ---- Build mapping table ----
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensembl_clean,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    SYMBOL     = SYMBOL
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

saveRDS(res_shrunk_list, "projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")














