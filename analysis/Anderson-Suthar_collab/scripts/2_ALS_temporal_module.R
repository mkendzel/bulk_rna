# ============================================================
# ALS Temporal Module score
# ============================================================
# ---- Library and Functions ----
library(tidyverse)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# ---- Load data ----

v <- load_checkpoint("voom_hSpS_Ctrx_batch", dir = "projects/Anderson-Suthar_collab/data/r_objects")
meta_filtered <- load_checkpoint("meta_filtered_hSpS_Ctrx", dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d30  <- load_checkpoint("res_d30_ALS_vs_Ctrl", dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d50  <- load_checkpoint("res_d50_ALS_vs_Ctrl", dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d75  <- load_checkpoint("res_d75_ALS_vs_Ctrl", dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d120 <- load_checkpoint("res_d120_ALS_vs_Ctrl", dir = "projects/Anderson-Suthar_collab/data/r_objects")

# ---- Determine genes/function for module score from d120 vs ctrl ----
# User parameters
fdr_cutoff  <- 0.05
lfc_cutoff  <- 1
res_input   <- res_d120

# Define module genes
module_genes <- res_input %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene") %>%
  filter(adj.P.Val < fdr_cutoff, abs(logFC) > lfc_cutoff)

# Summary 
cat("Total module genes:", nrow(module_genes), "\n")
cat("Upregulated in ALS:", sum(module_genes$logFC > 0), "\n")
cat("Downregulated in ALS:", sum(module_genes$logFC < 0), "\n")

# Top genes by absolute logFC
module_genes %>%
  arrange(desc(abs(logFC))) %>%
  dplyr::select(gene, logFC, AveExpr, adj.P.Val) %>%
  head(20) %>%
  print()

# Gene hits
als_genes <- data.frame(
  gene = c(
    "ANXA11", "ATXN1", "C9orf72", "CCNF", "CFAP410", "CHCHD10",
    "EPHA4", "FUS", "HFE", "HNRNPA1", "KIF5A", "NEK1", "NIPA1",
    "OPTN", "PFN1", "SCFD1", "SOD1", "TARDBP", "TBK1", "TUBA4A",
    "UBQLN2", "UNC13A", "VAPB", "VCP"
  ),
  als_evidence = c(
    "Definitive", "Strong", "Definitive", "Strong", "Strong", "Definitive",
    "Definitive", "Definitive", "Strong", "Definitive", "Definitive", "Definitive", "Strong",
    "Definitive", "Definitive", "Strong", "Definitive", "Definitive", "Definitive", "Strong",
    "Definitive", "Definitive", "Definitive", "Definitive"
  )
)

als_hits <- module_genes %>%
  inner_join(als_genes, by = "gene") %>%
  dplyr::select(gene, als_evidence, logFC, adj.P.Val) %>%
  arrange(desc(abs(logFC)))

cat("ALS genes in module:", nrow(als_hits), "of", nrow(module_genes), "total\n\n")
print(als_hits)


# Kegg als
# Fetch genes for ALS pathway directly from KEGG REST API
als_kegg_raw <- read.table(
  "https://rest.kegg.jp/link/hsa/hsa05014",
  sep    = "\t",
  header = FALSE,
  col.names = c("pathway", "gene")
)

# Strip "hsa:" prefix to get Entrez IDs
als_kegg_entrez <- gsub("hsa:", "", als_kegg_raw$gene)

# Map Entrez IDs to symbols
als_kegg_symbols <- mapIds(
  org.Hs.eg.db,
  keys      = als_kegg_entrez,
  column    = "SYMBOL",
  keytype   = "ENTREZID",
  multiVals = "first"
) %>%
  na.omit() %>%
  as.character()

cat("Genes in KEGG ALS pathway (hsa05014):", length(als_kegg_symbols), "\n")

# Overlap with module genes
kegg_hits <- module_genes %>%
  filter(gene %in% als_kegg_symbols)

cat("Module genes in KEGG ALS pathway:", nrow(kegg_hits), "\n")
print(kegg_hits %>% arrange(desc(abs(logFC))) %>% dplyr::select(gene, logFC, adj.P.Val))


# ---- Calculate module score ----
# Subset expression matrix to module genes
expr_module <- v$E[module_genes$gene, , drop = FALSE]

# Z-score each gene across samples
expr_scaled <- t(scale(t(expr_module)))

# Flip sign of downregulated genes so all contribute in the ALS-high direction
gene_sign <- sign(module_genes$logFC)
names(gene_sign) <- module_genes$gene
expr_directed <- expr_scaled * gene_sign[rownames(expr_scaled)]

# Average per sample
module_scores <- colMeans(expr_directed)

# Attach scores to metadata
score_df <- meta_filtered %>%
  mutate(module_score = module_scores[Admera_Health_ID])

# Summarize by Stage and Line
score_df %>%
  group_by(Stage, Line) %>%
  summarise(
    n          = n(),
    mean_score = mean(module_score),
    sd_score   = sd(module_score),
    .groups    = "drop"
  ) %>%
  print()

# ---- Module score for WNV ----
v_wnv <- v

module_genes_wnv <- module_genes %>%
  filter(gene %in% rownames(v_wnv$E))

cat("Module genes defined in ALS:", nrow(module_genes), "\n")
cat("Module genes present in WNV matrix:", nrow(module_genes_wnv), "\n")

# Subset WNV expression matrix to overlapping module genes
expr_module_wnv <- v_wnv$E[module_genes_wnv$gene, , drop = FALSE]

# Z-score each gene across WNV samples
expr_scaled_wnv <- t(scale(t(expr_module_wnv)))

# Orient each gene in the ALS-high direction using ALS-derived logFC sign
gene_sign <- sign(module_genes_wnv$logFC)
names(gene_sign) <- module_genes_wnv$gene
expr_directed_wnv <- expr_scaled_wnv * gene_sign[rownames(expr_scaled_wnv)]

# Average per sample
module_scores_wnv <- colMeans(expr_directed_wnv)

# Attach scores to WNV metadata
score_df_wnv <- meta_R2SDHF %>%
  mutate(module_score = module_scores_wnv[sample])

# Summarize by treatment and time
score_df_wnv %>%
  group_by(treatment, time) %>%
  summarise(
    n          = n(),
    mean_score = mean(module_score),
    sd_score   = sd(module_score),
    .groups    = "drop"
  ) %>%
  print()
