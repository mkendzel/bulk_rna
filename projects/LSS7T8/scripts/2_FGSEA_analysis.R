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

# =================================================
# ---- limma-voom GSVA and fgsea analyses ----
# =================================================

# ---- Limma FgSEA for LSS7T8 Liver ----
library(fgsea)
res_voom_liver_annotated <- load_checkpoint("res_voom_liver_annotated",
                dir = "projects/LSS7T8/data/r_objects")

# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

res_cl_vs_mock_liver   <- res_voom_liver_annotated[["treatmentcl"]]
res_dopc_vs_mock_liver <- res_voom_liver_annotated[["treatmentdopc"]]
res_dopc_vs_cl_liver   <- res_voom_liver_annotated[["treatment_dopc_vs_cl"]]

# Run for each treatment comparison
fgsea_cl_vs_mock_liver   <- run_gsea(res_cl_vs_mock_liver, pathways)
fgsea_dopc_vs_mock_liver <- run_gsea(res_dopc_vs_mock_liver, pathways)
fgsea_dopc_vs_cl_liver   <- run_gsea(res_dopc_vs_cl_liver, pathways)


save_checkpoint(fgsea_cl_vs_mock_liver,   "fgsea_cl_vs_mock_liver", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for cl vs mock in liver using voom limma results")

save_checkpoint(fgsea_dopc_vs_mock_liver, "fgsea_dopc_vs_mock_liver", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for dopc vs mock in liver using voom limma results")

save_checkpoint(fgsea_dopc_vs_cl_liver,   "fgsea_dopc_vs_cl_liver", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for dopc vs cl in liver using voom limma")


# ---- Limma FgSEA for LSS7T8 spleen ----
library(fgsea)

res_voom_spleen_annotated <- load_checkpoint("res_voom_spleen_annotated",
                dir = "projects/LSS7T8/data/r_objects")

# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

res_cl_vs_mock_spleen   <- res_voom_spleen_annotated[["treatmentcl"]]
res_dopc_vs_mock_spleen <- res_voom_spleen_annotated[["treatmentdopc"]]
res_dopc_vs_cl_spleen   <- res_voom_spleen_annotated[["treatment_dopc_vs_cl"]]

# Run for each treatment comparison
fgsea_cl_vs_mock_spleen   <- run_gsea(res_cl_vs_mock_spleen, pathways)
fgsea_dopc_vs_mock_spleen <- run_gsea(res_dopc_vs_mock_spleen, pathways)
fgsea_dopc_vs_cl_spleen   <- run_gsea(res_dopc_vs_cl_spleen, pathways)

save_checkpoint(fgsea_cl_vs_mock_spleen,   "fgsea_cl_vs_mock_spleen", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for cl vs mock in spleen using voom limma results")

save_checkpoint(fgsea_dopc_vs_mock_spleen, "fgsea_dopc_vs_mock_spleen", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for dopc vs mock in spleen using voom limma results")

save_checkpoint(fgsea_dopc_vs_cl_spleen,   "fgsea_dopc_vs_cl_spleen", dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA results for dopc vs cl in spleen using voom limma results")

# =================================================
# ---- DEseq 2 GSVA and fgsea analyses ----
# =================================================
# ---- GSVA for LSS7T8 Liver ----
# normalized counts with gene symbols as rownames
norm_counts_liver <- readRDS("projects/LSS7T8/data/r_objects/norm_counts_liver_symbols.rds")

# gsva parameters
gsva_par_liver <- gsvaParam(
  exprData  = norm_counts_liver,
  geneSets  = pathways,
  minSize   = 10,
  maxSize   = 500
)

#run GSVA
gsva_scores_liver <- gsva(gsva_par_liver, verbose = TRUE)

# Set up limma design and contrasts for GSVA scores
treatment_liver <- factor(colData(dds_liver)$treatment, levels = c("control", "cl", "dopc"))
design_liver <- model.matrix(~ 0 + treatment_liver)
colnames(design_liver) <- levels(treatment_liver)

#
contrasts_liver <- makeContrasts(
  cl_vs_control   = cl - control,
  dopc_vs_control = dopc - control,
  dopc_vs_cl      = dopc - cl,
  levels = design_liver
)

fit_liver <- lmFit(gsva_scores_liver, design_liver)
fit_liver <- contrasts.fit(fit_liver, contrasts_liver)

gssizes_liver <- rowSums(gsva_scores_liver != 0)
fit_liver <- eBayes(fit_liver, robust = TRUE, trend = gssizes_liver)

gsva_results_liver <- list(
  cl_vs_control   = topTable(fit_liver, coef = "cl_vs_control",   n = Inf, sort.by = "p"),
  dopc_vs_control = topTable(fit_liver, coef = "dopc_vs_control", n = Inf, sort.by = "p"),
  dopc_vs_cl      = topTable(fit_liver, coef = "dopc_vs_cl",      n = Inf, sort.by = "p")
)

saveRDS(gsva_scores_liver,  "projects/LSS7T8/data/r_objects/gsva_scores_liver.rds")
saveRDS(gsva_results_liver, "projects/LSS7T8/data/r_objects/gsva_results_liver.rds")


# ---- FgSEA for LSS7T8 Liver ----
res_shrunk_liver_annotated <- readRDS("projects/LSS7T8/data/r_objects/resShrink_liver_annotated.rds")
res_wald_liver_annotated   <- readRDS("projects/LSS7T8/data/r_objects/resWald_liver_annotated.rds")
# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)
res_wald_liver_annotated <- readRDS("projects/LSS7T8/data/r_objects/resWald_liver_annotated.rds")

# Helper to build ranked list and run fgsea
run_gsea <- function(tt, pathways, nPermSimple = 10000) {
  tt <- as.data.frame(tt)
  # Resolve duplicates by keeping the entry with highest baseMean
  n_before <- nrow(tt)
  tt <- tt[order(-tt$baseMean), ]
  tt <- tt[!duplicated(tt$symbol), ]
  n_dupes <- n_before - nrow(tt)
  message("Removed ", n_dupes, " duplicate gene names (kept highest baseMean)")

  ranks <- setNames(tt$stat, tt$symbol)
  ranks <- ranks[!is.na(names(ranks))]
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)

  res <- fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
  res <- res[order(res$pval), ]
  return(res)
}

res_cl_vs_mock_liver   <- res_wald_liver_annotated[["treatment_cl_vs_control"]]
res_dopc_vs_mock_liver <- res_wald_liver_annotated[["treatment_dopc_vs_control"]]
res_dopc_vs_cl_liver   <- res_wald_liver_annotated[["treatment_dopc_vs_cl"]]

fgsea_cl_vs_mock_liver   <- run_gsea(res_cl_vs_mock_liver, pathways)
fgsea_dopc_vs_mock_liver <- run_gsea(res_dopc_vs_mock_liver, pathways)
fgsea_dopc_vs_cl_liver   <- run_gsea(res_dopc_vs_cl_liver, pathways)

saveRDS(fgsea_cl_vs_mock_liver,   "projects/LSS7T8/data/r_objects/fgsea_cl_vs_mock_liver.rds")
saveRDS(fgsea_dopc_vs_mock_liver, "projects/LSS7T8/data/r_objects/fgsea_dopc_vs_mock_liver.rds")
saveRDS(fgsea_dopc_vs_cl_liver,   "projects/LSS7T8/data/r_objects/fgsea_dopc_vs_cl_liver.rds")

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

# DESEQ object creation
dds_spleen <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_spleen)),
  colData = info_spleen,
  design = ~ treatment
)

# (This object is only used to check if biomarkers were filtered out)
dds_spleen_unfiltered <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts_spleen)),
  colData = info_spleen,
  design = ~ treatment
)

#Subset low expression genes out using filterByExpr from edgeR
keep <- filterByExpr(dds_spleen, group = dds_spleen$treatment)
dds_spleen <- dds_spleen[keep, ]


## Check biomarkers ##
genes_of_interest <- c("Ifnb1", "Irf7", "Ifng", "Il6", "Cxcl10")

goi_ensembl <- na.omit(mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL %in% genes_of_interest])

counts_unfilt <- DESeq2::counts(dds_spleen_unfiltered)
counts_filt   <- DESeq2::counts(dds_spleen)
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
vsd_spleen <- DESeq2::vst(dds_spleen)
DESeq2::plotPCA(vsd_spleen, intgroup = "treatment")

# run DESeq2
dds_spleen <- DESeq2::DESeq(dds_spleen, minReplicatesForReplace = Inf)

#Grab coefficient names (except intercept)
coef_names_spleen <- DESeq2::resultsNames(dds_spleen)[-1]

# Wald test stat results (full comp)
res_wald_spleen <- setNames(
  lapply(coef_names_spleen, function(coef_name) {
    DESeq2::results(dds_spleen, name = coef_name, independentFiltering = FALSE)
  }),
  coef_names_spleen
)


# Add the specific contrast that DESeq2 doesn't automatically calculate (cl vs dopc)
res_wald_spleen[["treatment_dopc_vs_cl"]] <- DESeq2::results(
  dds_spleen,
  contrast = c("treatment", "dopc", "cl"),
  independentFiltering = FALSE
)

# Shrink coefficients using apeglm for all comparisons
res_shrunk_spleen <- setNames(
  lapply(coef_names_spleen, function(coef_name) {
    DESeq2::lfcShrink(dds_spleen, coef = coef_name, type = "apeglm")
  }),
  coef_names_spleen
)

# Shrink the specific contrast that DESeq2 doesn't automatically calculate (cl vs dopc) using ashr
res_shrunk_spleen[["treatment_dopc_vs_cl"]] <- DESeq2::lfcShrink(
  dds_spleen,
  res  = res_wald_spleen[["treatment_dopc_vs_cl"]],
  type = "ashr"
)

#Save/load DESeq object
saveRDS(dds_spleen,"projects/LSS7T8/data/r_objects/dds_spleen.rds")
saveRDS(DESeq2::counts(dds_spleen, normalized = TRUE), "projects/LSS7T8/data/r_objects/norm_counts_spleen.rds")
saveRDS(res_wald_spleen, "projects/LSS7T8/data/r_objects/res_wald_spleen.rds")
saveRDS(res_shrunk_spleen, "projects/LSS7T8/data/r_objects/res_shrunk_spleen.rds")

res_wald_spleen <- readRDS("projects/LSS7T8/data/r_objects/res_wald_spleen.rds")
res_shrunk_spleen <- readRDS("projects/LSS7T8/data/r_objects/res_shrunk_spleen.rds")
# ---- Annotated spleen----

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

res_shrunk_spleen_annotated <- purrr::map(res_shrunk_spleen, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_spleen, by = "ensembl_id") |>
    as_tibble()
})

res_wald_spleen_annotated <- purrr::map(res_wald_spleen, function(df) {
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map_spleen, by = "ensembl_id") |>
    as_tibble()
})

# Save Annotated
saveRDS(res_shrunk_spleen_annotated, "projects/LSS7T8/data/r_objects/resShrink_spleen_annotated.rds")
saveRDS(res_wald_spleen_annotated, "projects/LSS7T8/data/r_objects/resWald_spleen_annotated.rds")

res_shrunk_spleen_annotated <- readRDS("projects/LSS7T8/data/r_objects/resShrink_spleen_annotated.rds")
res_wald_spleen_annotated <- readRDS("projects/LSS7T8/data/r_objects/resWald_spleen_annotated.rds")

# ---- FgSEA for LSS7T8 Spleen ----


# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)
res_wald_spleen_annotated <- readRDS("projects/LSS7T8/data/r_objects/resWald_spleen_annotated.rds")

# Helper to build ranked list and run fgsea
run_gsea <- function(tt, pathways, nPermSimple = 10000) {
  tt <- as.data.frame(tt)
  # Resolve duplicates by keeping the entry with highest baseMean
  n_before <- nrow(tt)
  tt <- tt[order(-tt$baseMean), ]
  tt <- tt[!duplicated(tt$symbol), ]
  n_dupes <- n_before - nrow(tt)
  message("Removed ", n_dupes, " duplicate gene names (kept highest baseMean)")
  
  ranks <- setNames(tt$stat, tt$symbol)
  ranks <- ranks[!is.na(names(ranks))]
  ranks <- sort(ranks[!is.na(ranks)], decreasing = TRUE)
  
  res <- fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
  res <- res[order(res$pval), ]
  return(res)
}

res_cl_vs_mock_spleen   <- res_wald_spleen_annotated[["treatment_cl_vs_control"]]
res_dopc_vs_mock_spleen <- res_wald_spleen_annotated[["treatment_dopc_vs_control"]]
res_dopc_vs_cl_spleen   <- res_wald_spleen_annotated[["treatment_dopc_vs_cl"]]

fgsea_cl_vs_mock_spleen   <- run_gsea(res_cl_vs_mock_spleen, pathways)
fgsea_dopc_vs_mock_spleen <- run_gsea(res_dopc_vs_mock_spleen, pathways)
fgsea_dopc_vs_cl_spleen   <- run_gsea(res_dopc_vs_cl_spleen, pathways)

saveRDS(fgsea_cl_vs_mock_spleen,   "projects/LSS7T8/data/r_objects/fgsea_cl_vs_mock_spleen.rds")
saveRDS(fgsea_dopc_vs_mock_spleen, "projects/LSS7T8/data/r_objects/fgsea_dopc_vs_mock_spleen.rds")
saveRDS(fgsea_dopc_vs_cl_spleen,   "projects/LSS7T8/data/r_objects/fgsea_dopc_vs_cl_spleen.rds")

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