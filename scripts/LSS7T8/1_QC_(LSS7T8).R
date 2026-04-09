#### Joe's Data:Part 1 ####
# Code LSS7T8

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
library(TissueEnrich)
library(GSEABase)
library(edgeR)
# ----- Import data ----

# Geneset maps (load one based on model organism of interest)
mouse_gene_map <- readRDS("data/geneset/mouse_gene_map_ensembl_symbol.rds")
human_gene_map <- readRDS("data/geneset/human_gene_map_ensembl_symbol.rds")
#Import expression matrix in tsv format
expr <- read.delim(
  "data/LSS7T8/plasmidsaurus/LSS7T8-expression-matrix.tsv",
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

# remove low count samples (kept only as an example) THere are no low counts
# counts <- counts[, colnames(counts) != "cl_spleen_6"] # CS12
# counts <- counts[, colnames(counts) != "d_spleen_9"] # DS18

#Save/load checkpoint
saveRDS(counts, "data/LSS7T8/r_objects/plasmid_counts(all).rds")
counts <- readRDS("data/LSS7T8/r_objects/plasmid_counts(all).rds")

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

# ==================================
# ---- Tissue-specific DESeq 2 analyses ----
# ==================================
# ---- liver ----

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
res_wald_liver[["treatment_cl_vs_dopc"]] <- DESeq2::results(
  dds_liver,
  contrast = c("treatment", "cl", "dopc"),
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
res_shrunk_liver[["treatment_cl_vs_dopc"]] <- DESeq2::lfcShrink(
  dds_liver,
  res  = res_wald_liver[["treatment_cl_vs_dopc"]],
  type = "ashr"
)



#Save/load DESeq object
saveRDS(dds_liver,"data/LSS7T8/r_objects/dds_liver.rds")
saveRDS(DESeq2::counts(dds_liver, normalized = TRUE), "data/LSS7T8/r_objects/norm_counts_liver.rds")
saveRDS(res_wald_liver, "data/LSS7T8/r_objects/res_wald_liver.rds")
saveRDS(res_shrunk_liver, "data/LSS7T8/r_objects/res_shrunk_liver.rds")


# Annotated object
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

res_shrunk_liver <- purrr::map(res_shrunk_liver, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_liver, by = "ensembl_id") |>
    as_tibble()
})

res_wald_liver <- purrr::map(res_wald_liver, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_liver, by = "ensembl_id") |>
    as_tibble()
})

# Save Annotated
saveRDS(res_shrunk_liver, "data/LSS7T8/r_objects/plasmo_resShrink_liver_qcmin10_annotated.rds")
saveRDS(res_wald_liver, "data/LSS7T8/r_objects/plasmo_resWald_liver_qcmin10_annotated.rds")

res_shrunk_liver <- readRDS("data/LSS7T8/r_objects/plasmo_resShrink_liver_qcmin10_annotated.rds")
res_wald_liver <- readRDS("data/LSS7T8/r_objects/plasmo_resWald_liver_qcmin10_annotated.rds")

# Tissue check

vsd_liver <- assay(vst(dds_liver, blind = TRUE))

results_list_liver <- list()
for (s in colnames(vsd_liver)) {
  top_genes <- names(sort(vsd_liver[, s], decreasing = TRUE))[1:500]
  
  gs <- GeneSet(
    geneIds = top_genes,
    organism = "Mus Musculus",
    geneIdType = ENSEMBLIdentifier()
  )
  
  enrich <- teEnrichment(inputGenes = gs, rnaSeqDataset = 3)
  res <- as.data.frame(assay(enrich[[1]]))
  results_list_liver[[s]] <- res
}

for (s in names(results_list_liver)) {
  res <- results_list_liver[[s]]
  top_idx <- which.max(res$Tissue.Specific.Genes)
  cat(s, "->", rownames(res)[top_idx], 
      "(", res$Tissue.Specific.Genes[top_idx], "genes, p =", 
      10^(-res$Log10PValue[top_idx]), ")\n")
}


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
saveRDS(dds_spleen, "data/LSS7T8/r_objects/plasmo_dds_spleen_raw.rds")
dds_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_dds_spleen_raw.rds")

saveRDS(counts(dds_spleen, normalized = TRUE), "data/LSS7T8/r_objects/plasmo_counts_spleen_normalized.rds")
normalized_counts_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_counts_spleen_normalized.rds")

saveRDS(res_shrunk_spleen, "data/LSS7T8/r_objects/plasmo_resShrink_spleen_qcmin10.rds")
res_shrunk_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_resShrink_spleen_qcmin10.rds")

saveRDS(res_wald_spleen, "data/LSS7T8/r_objects/plasmo_resWald_spleen_qcmin10.rds")
res_wald_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_resWald_spleen_qcmin10.rds")

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
saveRDS(res_shrunk_spleen, "data/LSS7T8/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")
saveRDS(res_wald_spleen, "data/LSS7T8/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")

res_wald_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")
res_shrunk_spleen <- readRDS("data/LSS7T8/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")


# Tissue check

vsd_spleen <- assay(vst(dds_spleen, blind = TRUE))

# Average expression across all 9 samples
results_list_spleen <- list()

for (s in colnames(vsd_spleen)) {
  top_genes <- names(sort(vsd_spleen[, s], decreasing = TRUE))[1:500]
  
  gs <- GeneSet(
    geneIds = top_genes,
    organism = "Mus Musculus",
    geneIdType = ENSEMBLIdentifier()
  )
  
  enrich <- teEnrichment(inputGenes = gs, rnaSeqDataset = 3)
  res <- as.data.frame(assay(enrich[[1]]))
  results_list_spleen[[s]] <- res
}

for (s in names(results_list_spleen)) {
  res <- results_list_spleen[[s]]
  top_idx <- which.max(res$Tissue.Specific.Genes)
  cat(s, "->", rownames(res)[top_idx], 
      "(", res$Tissue.Specific.Genes[top_idx], "genes, p =", 
      10^(-res$Log10PValue[top_idx]), ")\n")
}
