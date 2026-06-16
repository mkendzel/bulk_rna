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
# ---- Tissue-specific DESeq 2 analyses ----
# ==========================================
# ------------------ liver -----------------

#Subset to liver
liver_samples <- sample_info$tissue == "liver"
counts_liver <- counts[, liver_samples]
info_liver <- sample_info[liver_samples, ]
info_liver$treatment <- relevel(droplevels(info_liver$treatment), ref = "control")

# DESEQ object creation
dds_liver <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_liver)),
  colData = info_liver,
  design = ~ treatment
)

# (This object is only used to check if biomarkers were filtered out)
dds_liver_unfiltered <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_liver)),
  colData = info_liver,
  design = ~ treatment
)

#Subset low expression genes out using filterByExpr from edgeR
keep <- filterByExpr(dds_liver, group = dds_liver$treatment)
dds_liver <- dds_liver[keep, ]


## Check biomarkers ##
genes_of_interest <- c("Ifnb1", "Irf7", "Ifng", "Il6", "Cxcl10")

goi_ensembl <- na.omit(mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL %in% genes_of_interest])

counts_unfilt <- DESeq2::counts(dds_liver_unfiltered)
counts_filt   <- DESeq2::counts(dds_liver)
cpm_unfilt    <- edgeR::cpm(counts_unfilt)
cpm_filt      <- edgeR::cpm(counts_filt)

goi_present_unfilt <- goi_ensembl[goi_ensembl %in% rownames(counts_unfilt)]
goi_present_filt   <- goi_ensembl[goi_ensembl %in% rownames(counts_filt)]

summarise_goi <- function(ensembl_ids, counts_mat, cpm_mat, gene_map) {
  data.frame(
    symbol        = gene_map$SYMBOL[match(ensembl_ids, gene_map$ENSEMBL)],
    mean_count    = rowMeans(counts_mat[ensembl_ids, , drop = FALSE]),
    median_count  = apply(counts_mat[ensembl_ids, , drop = FALSE], 1, median),
    max_count     = apply(counts_mat[ensembl_ids, , drop = FALSE], 1, max),
    mean_cpm      = rowMeans(cpm_mat[ensembl_ids, , drop = FALSE]),
    median_cpm    = apply(cpm_mat[ensembl_ids, , drop = FALSE], 1, median),
    n_samples_gt0 = rowSums(counts_mat[ensembl_ids, , drop = FALSE] > 0)
  )
}

goi_summary_unfilt <- summarise_goi(goi_present_unfilt, counts_unfilt, cpm_unfilt, mouse_gene_map)
goi_summary_filt   <- summarise_goi(goi_present_filt,   counts_filt,   cpm_filt,   mouse_gene_map)

lost <- goi_present_unfilt[!goi_present_unfilt %in% goi_present_filt]

list(
  before   = goi_summary_unfilt,
  after    = goi_summary_filt,
  filtered = if (length(lost) > 0) mouse_gene_map$SYMBOL[match(lost, mouse_gene_map$ENSEMBL)] else "all genes survived"
)
##

# variance stabilization for PCA
vsd_liver <- DESeq2::vst(dds_liver)
DESeq2::plotPCA(vsd_liver, intgroup = "treatment")

# run DESeq2
dds_liver <- DESeq2::DESeq(dds_liver, minReplicatesForReplace = Inf)

#Grab coefficient names (except intercept)
coef_names_liver <- DESeq2::resultsNames(dds_liver)[-1]

# Wald test stat results (full comp)
res_wald_liver <- setNames(
  lapply(coef_names_liver, function(coef_name) {
    DESeq2::results(dds_liver, name = coef_name, independentFiltering = FALSE)
  }),
  coef_names_liver
)


# Add the specific contrast that DESeq2 doesn't automatically calculate (cl vs dopc)
res_wald_liver[["treatment_dopc_vs_cl"]] <- DESeq2::results(
  dds_liver,
  contrast = c("treatment", "dopc", "cl"),
  independentFiltering = FALSE
)

# Shrink coefficients using apeglm for all comparisons
res_shrunk_liver <- setNames(
  lapply(coef_names_liver, function(coef_name) {
    DESeq2::lfcShrink(dds_liver, coef = coef_name, type = "apeglm")
  }),
  coef_names_liver
)

# Shrink the specific contrast that DESeq2 doesn't automatically calculate (cl vs dopc) using ashr
res_shrunk_liver[["treatment_dopc_vs_cl"]] <- DESeq2::lfcShrink(
  dds_liver,
  res  = res_wald_liver[["treatment_dopc_vs_cl"]],
  type = "ashr"
)

#Save/load DESeq object
saveRDS(dds_liver,"projects/LSS7T8/data/r_objects/dds_liver.rds")
saveRDS(DESeq2::counts(dds_liver, normalized = TRUE), "projects/LSS7T8/data/r_objects/norm_counts_liver.rds")
saveRDS(res_wald_liver, "projects/LSS7T8/data/r_objects/res_wald_liver.rds")
saveRDS(res_shrunk_liver, "projects/LSS7T8/data/r_objects/res_shrunk_liver.rds")

# -------------- Annotated liver------------
# Load objects
dds_liver        <- readRDS("projects/LSS7T8/data/r_objects/dds_liver.rds")
res_shrunk_liver <- readRDS("projects/LSS7T8/data/r_objects/res_shrunk_liver.rds")
res_wald_liver   <- readRDS("projects/LSS7T8/data/r_objects/res_wald_liver.rds")

all_ensembl_liver <- res_shrunk_liver |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

gene_map_liver <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys    = all_ensembl_liver,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)

gene_map_liver <- gene_map_liver |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    symbol     = SYMBOL
  )

res_shrunk_liver_annotated <- purrr::map(res_shrunk_liver, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_liver, by = "ensembl_id") |>
    as_tibble()
})

res_wald_liver_annotated <- purrr::map(res_wald_liver, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_liver, by = "ensembl_id") |>
    as_tibble()
})

# Build symbol-annotated normalized count matrix
norm_counts_liver <- counts(dds_liver, normalized = TRUE)
rownames(norm_counts_liver) <- sub("\\..*$", "", rownames(norm_counts_liver))

anno <- gene_map_liver |>
  filter(!is.na(symbol), !duplicated(ensembl_id))

norm_counts_liver <- norm_counts_liver[rownames(norm_counts_liver) %in% anno$ensembl_id, ]
symbol_map <- setNames(anno$symbol, anno$ensembl_id)
rownames(norm_counts_liver) <- symbol_map[rownames(norm_counts_liver)]
norm_counts_liver <- norm_counts_liver[!is.na(rownames(norm_counts_liver)), ]

# For duplicate symbols, keep the row with highest mean expression
norm_counts_liver <- norm_counts_liver[order(-rowMeans(norm_counts_liver)), ]
norm_counts_liver <- norm_counts_liver[!duplicated(rownames(norm_counts_liver)), ]

# Save all
saveRDS(res_shrunk_liver_annotated, "projects/LSS7T8/data/r_objects/resShrink_liver_annotated.rds")
saveRDS(res_wald_liver_annotated,   "projects/LSS7T8/data/r_objects/resWald_liver_annotated.rds")
saveRDS(norm_counts_liver,          "projects/LSS7T8/data/r_objects/norm_counts_liver_symbols.rds")







