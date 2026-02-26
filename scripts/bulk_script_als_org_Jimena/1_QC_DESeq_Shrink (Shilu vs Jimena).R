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
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

##### Jimina's data
#### Expression matrix
# Define base path
base_path <- "data/Anderson-Suthar_collab/Als_org_Jimena"

# Import expression matrix
expr <- read_excel(
  file.path(base_path, "102025_lines1-6,c9_vs_ctrl,d30,50,75,120_allfeature_counts.xlsx")
)

# rowname set up
expr <- as.data.frame(expr)
rownames(expr) <- expr$Geneid
expr <- expr[, -1]

# colname formating
colnames(expr) <- sub("^\\./", "", colnames(expr))
colnames(expr) <- sub("_Aligned.*$", "", colnames(expr))


#### Metadata
# Import metadata
metadata <- read_excel(
  file.path(base_path, "BulkSeqTargetList_tosend.xlsx")
)

metadata <- as.data.frame(metadata)

# ---- Set up and subset datasets ----

# Filter metadata for hSPS organoids with Ctrx coating
meta_filtered <- metadata[metadata$Organoid == "hSPS" & metadata$Coating == "Ctrx"]

# Keep only relevant columns
meta_filtered <- meta_filtered[, c("Admera_Helath_ID", "Patient", "AgeOrg")]

# Subset expression matrix to only matching samples
expr_filtered <- expr[, colnames(expr) %in% meta_filtered$Admera_Helath_ID]

# Make sure they're in the same order (critical for DESeq2)
meta_filtered <- meta_filtered[match(colnames(expr_filtered), meta_filtered$Admera_Helath_ID), ]

# Confirm order matches
all(colnames(expr_filtered) == meta_filtered$Admera_Helath_ID)

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
counts <- readRDS("data/R2SDHF/r_objects/plasmo_counts_noMAD.rds")

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
# saveRDS(dds, "data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")
dds <- readRDS("data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#Save/load Shrunk Data
saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo_resShrink_noMAD_qcmin10.rds")
res_shrunk_list <- readRDS("data/R2SDHF/Plasmo_resShrink_noMAD_qcmin10.rds")


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

saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("data/R2SDHF/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")














