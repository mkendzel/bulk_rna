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
## Primary import: expression matrix (TSV)
# EDIT: point this at your project's expression matrix
expr <- read.delim(
  "projects/PROJECT_NAME/data/PROJECT_NAME-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## Alternative import: per-contrast DE tables instead of a count matrix
# de_files <- list.files("projects/PROJECT_NAME/data", pattern = "^DE_", full.names = TRUE)
#
# de_list <- setNames(
#   lapply(de_files, function(f) {
#     if (grepl("\\.csv$", f)) read.csv(f) else read_excel(f)
#   }),
#   gsub("\\.(xlsx|csv)$", "", basename(de_files))
# )
#
# names(de_list)



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
  "SAMPLE_1"  = "Control_0_1",
  "SAMPLE_2"  = "Control_0_2",
  "SAMPLE_3"  = "Control_0_3",
  "SAMPLE_4"  = "TreatmentA_1_1",
  "SAMPLE_5"  = "TreatmentA_1_2",
 # etc.
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

# Drop low-count samples based on qc_summary/qc_list above
counts <- counts[, colnames(counts) != "TreatmentB_1_1"]
counts <- counts[, colnames(counts) != "TreatmentB_1_2"]


# Save/load checkpoint
# save_checkpoint(counts, "counts_filtered", dir = "projects/YOUR_PROJECT_NAME/data/r_objects")
counts <- load_checkpoint("counts_filtered", dir = "projects/YOUR_PROJECT_NAME/data/r_objects")
 
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
  colData   = colData,
  design    = ~ condition
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
coef_names_vs_ref <- coef_names[grep(paste0("vs_", control_condition), coef_names)]

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
dds <- load_checkpoint("dds_filtered", dir = "projects/PROJECT_NAME/data/r_objects")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

# Save/load normalized counts
save_checkpoint(normalized_counts, "normalized_counts", dir = "projects/PROJECT_NAME/data/r_objects")
normalized_counts <- load_checkpoint("normalized_counts", dir = "projects/PROJECT_NAME/data/r_objects")

# Save/load Shrunk Data
save_checkpoint(res_shrunk_list, "res_shrunk_list", dir = "projects/PROJECT_NAME/data/r_objects")
res_shrunk_list <- load_checkpoint("res_shrunk_list", dir = "projects/PROJECT_NAME/data/r_objects")


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

save_checkpoint(res_shrunk_list, "res_shrunk_list_annotated", dir = "projects/PROJECT_NAME/data/r_objects")
res_shrunk_list <- load_checkpoint("res_shrunk_list_annotated", dir = "projects/PROJECT_NAME/data/r_objects")