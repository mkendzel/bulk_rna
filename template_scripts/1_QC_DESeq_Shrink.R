# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)
library(readxl)
# ----- Import data ----
## Primary import: expression matrix (TSV)
# EDIT: point this at your project's expression matrix
expr <- read.delim(
  "projects/PROJECT_NAME/data/PROJECT_NAME-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## Alternative import: pre-computed per-contrast DE tables
# Use this instead of (or alongside) the matrix above if your data arrives as
# one DE results file per comparison (xlsx/csv) rather than a raw count matrix.
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


# rename columns
# EDIT: this is the key project-specific step. LHS = the raw column names in
# YOUR expression matrix (run `colnames(counts)` first to see them — e.g. a
# Plasmidsaurus matrix has columns like CODE_1, CODE_2, ...). RHS = meaningful
# `condition_timepoint_replicate` names that you choose to describe your design.
# The example below shows a control + 3-treatment, 3-timepoint, 3-replicate
# design (one of the more complex shapes you'll encounter) — trim/extend rows
# and rename groups to match your actual experiment.
rename_map <- c(
  "SAMPLE_1"  = "Control_0_1",
  "SAMPLE_2"  = "Control_0_2",
  "SAMPLE_3"  = "Control_0_3",
  "SAMPLE_4"  = "TreatmentA_1_1",
  "SAMPLE_5"  = "TreatmentA_1_2",
  "SAMPLE_6"  = "TreatmentA_1_3",
  "SAMPLE_7"  = "TreatmentC_1_1",
  "SAMPLE_8"  = "TreatmentC_1_2",
  "SAMPLE_9"  = "TreatmentC_1_3",
  "SAMPLE_10" = "TreatmentB_1_1",
  "SAMPLE_11" = "TreatmentB_1_2",
  "SAMPLE_12" = "TreatmentB_1_3",
  "SAMPLE_13" = "TreatmentA_2_1",
  "SAMPLE_14" = "TreatmentA_2_2",
  "SAMPLE_15" = "TreatmentA_2_3",
  "SAMPLE_16" = "TreatmentC_2_1",
  "SAMPLE_17" = "TreatmentC_2_2",
  "SAMPLE_18" = "TreatmentC_2_3",
  "SAMPLE_19" = "TreatmentB_2_1",
  "SAMPLE_20" = "TreatmentB_2_2",
  "SAMPLE_21" = "TreatmentB_2_3",
  "SAMPLE_22" = "TreatmentA_3_1",
  "SAMPLE_23" = "TreatmentA_3_2",
  "SAMPLE_24" = "TreatmentA_3_3",
  "SAMPLE_25" = "TreatmentC_3_1",
  "SAMPLE_26" = "TreatmentC_3_2",
  "SAMPLE_27" = "TreatmentC_3_3",
  "SAMPLE_28" = "TreatmentB_3_1",
  "SAMPLE_29" = "TreatmentB_3_2",
  "SAMPLE_30" = "TreatmentB_3_3"
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
# EDIT: this is a decision point — look at qc_summary / qc_list above (total
# reads, detected genes, percent of library) and decide which samples (if any)
# should be dropped before fitting the model. The example below shows an entire
# treatment group (TreatmentB) failing QC and being removed at every timepoint —
# this is exactly the kind of mid-pipeline decision the QC step is meant to
# surface. Replace these names with whatever samples your own QC flags; don't
# copy them verbatim.
counts <- counts[, colnames(counts) != "TreatmentB_1_1"]
counts <- counts[, colnames(counts) != "TreatmentB_1_2"]
counts <- counts[, colnames(counts) != "TreatmentB_1_3"]
counts <- counts[, colnames(counts) != "TreatmentB_3_1"]
counts <- counts[, colnames(counts) != "TreatmentB_3_2"]
counts <- counts[, colnames(counts) != "TreatmentB_3_3"]
counts <- counts[, colnames(counts) != "TreatmentB_2_1"]
counts <- counts[, colnames(counts) != "TreatmentB_2_2"]
counts <- counts[, colnames(counts) != "TreatmentB_2_3"]

#Save/load checkpoint
# Filename encodes the QC decision above (which samples were dropped) so later
# scripts/analyses stay traceable — adjust the tag to reflect your own decision.
# saveRDS(counts, "projects/PROJECT_NAME/data/r_objects/counts_filtered.rds")
counts <- readRDS("projects/PROJECT_NAME/data/r_objects/counts_filtered.rds")

# ---- set up deseq2 object ----
sample_names <- colnames(counts)

# Derive `condition` directly from the `condition_replicate` names chosen in
# rename_map above, by stripping the trailing replicate number. This works for
# any treatment/timepoint scheme as long as the convention is
# `<condition>_<replicate>` with a numeric replicate suffix — adjust the regex
# if your replicate isn't a trailing integer (e.g. "_rep1").
colData <- data.frame(sample = sample_names) |>
  dplyr::mutate(
    condition = sub("_[0-9]+$", "", sample)
  )

rownames(colData) <- colData$sample

# EDIT: set this to your control/reference condition so all contrasts are
# computed against it.
control_condition <- "Control_0"
colData$condition <- relevel(factor(colData$condition), ref = control_condition)

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = colData,
  design    = ~ condition
)

# Remove low-information genes to reduce noise and speed up fitting
# EDIT: these are typical starting points (>= 10 counts in at least N samples,
# where N is often the size of your smallest group) — revisit based on your
# own sample sizes and sequencing depth.
keep <- rowSums(DESeq2::counts(dds) >= 10) >= 5
dds <- dds[keep, ]

# Transform counts for sample-level exploration (PCA)
vsd <- DESeq2::vst(dds)
DESeq2::plotPCA(vsd, intgroup = "condition")

# Fit the DESeq2 model; disable outlier replacement for consistent behavior across conditions
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Inspect the coefficient names to confirm which contrasts exist
coef_names <- DESeq2::resultsNames(dds)

# Keep only coefficients that represent contrasts against the control reference
coef_names_vs_ref <- coef_names[grep(paste0("vs_", control_condition), coef_names)]

# Shrink log2 fold-changes with apeglm for each vs-control coefficient
res_shrunk_list <- setNames(
  lapply(coef_names_vs_ref, function(coef_name) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  }),
  coef_names_vs_ref
)

# Quick check of unshrunk MA plot for the default results
res <- DESeq2::results(dds)
DESeq2::plotMA(res, ylim = c(-5, 5))


#Save/load DESeq object
# saveRDS(dds, "projects/PROJECT_NAME/data/r_objects/dds_filtered.rds")
dds <- readRDS("projects/PROJECT_NAME/data/r_objects/dds_filtered.rds")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#Save/load Shrunk Data
saveRDS(res_shrunk_list, "projects/PROJECT_NAME/results/resShrink.rds")
res_shrunk_list <- readRDS("projects/PROJECT_NAME/results/resShrink.rds")


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
# Default below is human (org.Hs.eg.db). For mouse, swap for org.Mm.eg.db
# (the genesets/ folder also ships a pre-built mouse_gene_map_ensembl_symbol.rds
# you can read directly instead of querying AnnotationDbi).
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

saveRDS(res_shrunk_list, "projects/PROJECT_NAME/results/resShrink_annotated.rds")
res_shrunk_list <- readRDS("projects/PROJECT_NAME/results/resShrink_annotated.rds")














