#### Joe's Data:Part 1 ####
# Code LSS7T8

# ---- Libraries ----
library(dplyr)
library(DESeq2)
library(org.Mm.eg.db)
library(tibble)
library(TissueEnrich)
library(GSEABase)
library(edgeR)
library(fgsea)
library(GSVA)
library(limma)
# ----- Import data ----

#import all helper functions from folder R/
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# Geneset maps (load one based on model organism of interest)
mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")
# human_gene_map <- readRDS("genesets/human_gene_map_ensembl_symbol.rds")
#Import expression matrix in tsv format
expr <- read.delim(
  "projects/LSS7T8/data/plasmidsaurus/LSS7T8-expression-matrix.tsv",
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
  "LSS7T8_1"  = "control_liver_1",
  "LSS7T8_2"  = "control_liver_2",
  "LSS7T8_3"  = "control_liver_3",
  "LSS7T8_4"  = "control_spleen_1",
  "LSS7T8_5"  = "control_spleen_2",
  "LSS7T8_6"  = "control_spleen_3",
  "LSS7T8_7"  = "cl_liver_4",
  "LSS7T8_8"  = "cl_liver_5",
  "LSS7T8_9"  = "cl_liver_6",
  "LSS7T8_10" = "cl_spleen_4",
  "LSS7T8_11" = "cl_spleen_5",
  "LSS7T8_12" = "cl_spleen_6",
  "LSS7T8_13" = "dopc_liver_7",
  "LSS7T8_14" = "dopc_liver_8",
  "LSS7T8_15" = "dopc_liver_9",
  "LSS7T8_16" = "dopc_spleen_7",
  "LSS7T8_17" = "dopc_spleen_8",
  "LSS7T8_18" = "dopc_spleen_9"
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

# Notes from QC:
# There were no samples requiring removal in this data set #

#Save/load checkpoint
save_checkpoint(counts, "plasmid_counts", notes = "no samples removed after QC",
                dir = "projects/LSS7T8/data/r_objects")

list_checkpoints("plasmid_counts", dir = "projects/LSS7T8/data/r_objects")
counts <- load_checkpoint("plasmid_counts", dir = "projects/LSS7T8/data/r_objects")

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

# Save sample info
save_checkpoint(sample_info, "sample_info", notes = "parsed from column names",
                dir = "projects/LSS7T8/data/r_objects", lines = c(115:126))


# =============================================
# ---- Tissue-specific Limma-voom analyses ----
# =============================================
# ----------------- Liver ---------------------
counts <- load_checkpoint("plasmid_counts", dir = "projects/LSS7T8/data/r_objects")
sample_info <- load_checkpoint("sample_info", dir = "projects/LSS7T8/data/r_objects")
mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")

# Subset to liver
liver_samples <- sample_info$tissue == "liver"
counts_liver <- counts[, liver_samples]
info_liver <- sample_info[liver_samples, ]
info_liver$treatment <- relevel(droplevels(info_liver$treatment), ref = "control")

# Create DGEList and filter
dge <- DGEList(counts = round(as.matrix(counts_liver)), group = info_liver$treatment)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Design matrix and voom
design <- model.matrix(~ treatment, data = info_liver)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract results for reference-level comparisons
coef_names <- colnames(design)[-1]
res_voom_liver <- setNames(
  lapply(coef_names, function(coef) {
    topTable(fit, coef = coef, number = Inf, sort.by = "none")
  }),
  coef_names
)

# Add DOPC vs CL contrast
contr_matrix <- makeContrasts(
  treatmentdopc_vs_cl = treatmentdopc - treatmentcl,
  levels = design
)
fit2 <- contrasts.fit(fit, contr_matrix)
fit2 <- eBayes(fit2)
res_voom_liver[["treatment_dopc_vs_cl"]] <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")

# Annotate with gene symbols
res_voom_liver_annotated <- lapply(res_voom_liver, function(df) {
  df %>%
    rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) %>%
    left_join(
      mouse_gene_map %>% distinct(ENSEMBL, .keep_all = TRUE),
      by = c("ensembl_id" = "ENSEMBL")
    ) %>%
    as_tibble()
})

# Save
save_checkpoint(v, "voom_liver", notes = "voom object for liver samples", 
                dir = "projects/LSS7T8/data/r_objects", lines = c(147:155))

save_checkpoint(res_voom_liver, "res_voom_liver",
                notes = "limma-voom DE results for liver; named list of topTables (each treatment vs control + dopc-vs-cl contrast)",
                dir = "projects/LSS7T8/data/r_objects", lines = c(159:175))

save_checkpoint(res_voom_liver_annotated, "res_voom_liver_annotated",
                notes = "Gene Symbol annotated limma-voom DE results for liver; named list of topTables (each treatment vs control + dopc-vs-cl contrast)",
                dir = "projects/LSS7T8/data/r_objects", lines = c(177:187))

# ------------------- Spleen ---------------------
counts <- load_checkpoint("plasmid_counts", dir = "projects/LSS7T8/data/r_objects")
sample_info <- load_checkpoint("sample_info", dir = "projects/LSS7T8/data/r_objects")
mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")

# Subset to spleen
spleen_samples <- sample_info$tissue == "spleen"
counts_spleen <- counts[, spleen_samples]
info_spleen <- sample_info[spleen_samples, ]
info_spleen$treatment <- relevel(droplevels(info_spleen$treatment), ref = "control")

# Create DGEList and filter
dge <- DGEList(counts = round(as.matrix(counts_spleen)), group = info_spleen$treatment)
keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

# Design matrix and voom
design <- model.matrix(~ treatment, data = info_spleen)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# Extract results for reference-level comparisons
coef_names <- colnames(design)[-1]
res_voom_spleen <- setNames(
  lapply(coef_names, function(coef) {
    topTable(fit, coef = coef, number = Inf, sort.by = "none")
  }),
  coef_names
)

# Add CL vs DOPC contrast
contr_matrix <- makeContrasts(
  treatmentdopc_vs_cl = treatmentdopc - treatmentcl,
  levels = design
)
fit2 <- contrasts.fit(fit, contr_matrix)
fit2 <- eBayes(fit2)
res_voom_spleen[["treatment_dopc_vs_cl"]] <- topTable(fit2, coef = 1, number = Inf, sort.by = "none")

# Annotate with gene symbols
res_voom_spleen_annotated <- lapply(res_voom_spleen, function(df) {
  df %>%
    rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) %>%
    left_join(
      mouse_gene_map %>% distinct(ENSEMBL, .keep_all = TRUE),
      by = c("ensembl_id" = "ENSEMBL")
    ) %>%
    as_tibble()
})

# Save
save_checkpoint(v, "voom_spleen", notes = "voom object for spleen samples", 
                dir = "projects/LSS7T8/data/r_objects", lines = c(206:220))

save_checkpoint(res_voom_spleen, "res_voom_spleen",
                notes = "limma-voom DE results for spleen; named list of topTables (each treatment vs control + dopc-vs-cl contrast)",
                dir = "projects/LSS7T8/data/r_objects", lines = c(221:240))

save_checkpoint(res_voom_spleen_annotated, "res_voom_spleen_annotated",
                notes = "Gene Symbol annotated limma-voom DE results for spleen; named list of topTables (each treatment vs control + dopc-vs-cl contrast)",
                dir = "projects/LSS7T8/data/r_objects", lines = c(242:252))

# ==========================================
