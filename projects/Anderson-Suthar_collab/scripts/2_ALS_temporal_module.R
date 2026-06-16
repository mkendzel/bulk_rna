# ============================================================
# ALS Temporal Module: Genes Whose Disease Signature Strengthens
# ============================================================
# Starting from the per-timepoint DE results (ALS vs control) already
# computed in script 1, this script identifies genes where the ALS-vs-control
# logFC increases monotonically from d30 → d50 → d75 → d120.

library(tidyverse)
library(ggplot2)
library(pheatmap)
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)

source("R/checkpoint.R")
source("R/checkpoint_meta.R")

# ---- 1. Load saved DE results ----

res_d30  <- readRDS("projects/Anderson-Suthar_collab/data/r_objects/res_d30_ALS_vs_Ctrl(Batch).rds")
res_d50  <- readRDS("projects/Anderson-Suthar_collab/data/r_objects/res_d50_ALS_vs_Ctrl(Batch).rds")
res_d75  <- readRDS("projects/Anderson-Suthar_collab/data/r_objects/res_d75_ALS_vs_Ctrl(Batch).rds")
res_d120 <- readRDS("projects/Anderson-Suthar_collab/data/r_objects/res_d120_ALS_vs_Ctrl(Batch).rds")

# All four objects must share the same gene universe
stopifnot(identical(rownames(res_d30), rownames(res_d50)),
          identical(rownames(res_d30), rownames(res_d75)),
          identical(rownames(res_d30), rownames(res_d120)))

genes <- rownames(res_d30)
message(length(genes), " genes in the background set")

# ---- 2. Build logFC-across-time matrix ----

lfc_mat <- data.frame(
  gene_id    = genes,
  lfc_d30    = res_d30[genes,  "logFC"],
  lfc_d50    = res_d50[genes,  "logFC"],
  lfc_d75    = res_d75[genes,  "logFC"],
  lfc_d120   = res_d120[genes, "logFC"],
  adj_p_d30  = res_d30[genes,  "adj.P.Val"],
  adj_p_d50  = res_d50[genes,  "adj.P.Val"],
  adj_p_d75  = res_d75[genes,  "adj.P.Val"],
  adj_p_d120 = res_d120[genes, "adj.P.Val"],
  row.names  = genes,
  stringsAsFactors = FALSE
)

lfc_mat$minFDR <- pmin(
  lfc_mat$adj_p_d30, lfc_mat$adj_p_d50,
  lfc_mat$adj_p_d75, lfc_mat$adj_p_d120
)

# Spearman rho of logFC vs time — useful as a continuous ranking metric
timepoints <- c(30, 50, 75, 120)
lfc_mat$spearman_rho <- apply(
  lfc_mat[, c("lfc_d30", "lfc_d50", "lfc_d75", "lfc_d120")], 1,
  function(x) cor(x, timepoints, method = "spearman")
)

# ---- 3. Filter: strictly monotonic upregulation ----
# Criteria:
#   logFC increases at every step (d30 < d50 < d75 < d120)
#   Net positive at the endpoint (d120 > 0)
#   Significant in at least one stage (minFDR < 0.05)

als_up_module <- lfc_mat %>%
  filter(
    lfc_d50  > lfc_d30,
    lfc_d75  > lfc_d50,
    lfc_d120 > lfc_d75,
    lfc_d120 > 0,
    minFDR   < 0.05
  ) %>%
  arrange(desc(lfc_d120))

message(nrow(als_up_module), " genes pass strict monotonic filter")

# ---- 4. Gene symbol mapping ----
# If gene IDs look like Ensembl IDs, map to HGNC symbols.
# The existing pipeline uses rownames directly for GSEA (symbol-based pathways),
# so if no mapping is needed the gene_id column already contains symbols.

if (grepl("^ENSG", als_up_module$gene_id[1])) {
  message("Ensembl IDs detected — mapping to HGNC symbols via org.Hs.eg.db")

  clean_ids <- sub("\\..*$", "", als_up_module$gene_id)  # strip version suffix
  sym_map <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys    = clean_ids,
    keytype = "ENSEMBL",
    columns = c("SYMBOL")
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE)

  als_up_module <- als_up_module %>%
    mutate(ensembl_clean = sub("\\..*$", "", gene_id)) %>%
    left_join(sym_map, by = c("ensembl_clean" = "ENSEMBL")) %>%
    rename(symbol = SYMBOL) %>%
    filter(!is.na(symbol)) %>%
    arrange(desc(lfc_d120)) %>%
    distinct(symbol, .keep_all = TRUE)

} else {
  message("Gene symbols detected — no ID mapping needed")
  als_up_module$symbol <- als_up_module$gene_id
}

message(nrow(als_up_module), " genes after symbol resolution")

# ---- 5. Save the module ----

dir.create("projects/Anderson-Suthar_collab/results", recursive = TRUE, showWarnings = FALSE)

save_checkpoint(
  als_up_module,
  name = "als_up_module",
  dir  = "projects/Anderson-Suthar_collab/data/r_objects"
)

write.csv(
  als_up_module,
  "projects/Anderson-Suthar_collab/results/als_up_module_genes.csv",
  row.names = FALSE
)

# ---- 6. Visualizations ----

stage_colors <- c(d30 = "#C6DBEF", d50 = "#6BAED6", d75 = "#2171B5", d120 = "#08306B")

# A. Line plot: mean logFC trajectory of the module
long_lfc <- als_up_module %>%
  dplyr::select(gene_id, lfc_d30, lfc_d50, lfc_d75, lfc_d120) %>%
  pivot_longer(
    cols      = starts_with("lfc_"),
    names_to  = "stage",
    values_to = "logFC"
  ) %>%
  mutate(stage = factor(
    recode(stage, lfc_d30 = "d30", lfc_d50 = "d50", lfc_d75 = "d75", lfc_d120 = "d120"),
    levels = c("d30", "d50", "d75", "d120")
  ))

module_summary <- long_lfc %>%
  group_by(stage) %>%
  summarise(
    mean_lfc = mean(logFC),
    se_lfc   = sd(logFC) / sqrt(n()),
    .groups  = "drop"
  )

p_line <- ggplot(module_summary, aes(x = stage, y = mean_lfc, group = 1)) +
  geom_ribbon(aes(ymin = mean_lfc - se_lfc, ymax = mean_lfc + se_lfc),
              fill = "#CB181D", alpha = 0.15) +
  geom_line(color = "#CB181D", linewidth = 1.2) +
  geom_point(aes(fill = stage), shape = 21, size = 4, color = "white",
             fill = "#CB181D") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x     = "Organoid Stage",
    y     = "Mean logFC (ALS vs Control) ± SE",
    title = paste0("ALS Upward Module (n = ", nrow(als_up_module), " genes)"),
    subtitle = "Disease signature strengthens through organoid development"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")

ggsave(
  "projects/Anderson-Suthar_collab/figures/als_up_module_trajectory.pdf",
  p_line, width = 5, height = 4
)

# B. Heatmap: module logFC across stages (top 80 by d120 logFC)
n_heatmap <- min(nrow(als_up_module), 80)
top_genes  <- als_up_module %>% slice_max(order_by = lfc_d120, n = n_heatmap) %>% pull(symbol)

lfc_heatmap <- as.matrix(
  als_up_module[als_up_module$symbol %in% top_genes,
                c("lfc_d30", "lfc_d50", "lfc_d75", "lfc_d120")]
)
rownames(lfc_heatmap) <- als_up_module$symbol[als_up_module$symbol %in% top_genes]
colnames(lfc_heatmap) <- c("d30", "d50", "d75", "d120")

abs_max <- max(abs(lfc_heatmap))

pheatmap(
  lfc_heatmap,
  cluster_cols  = FALSE,
  cluster_rows  = TRUE,
  scale         = "none",
  color         = colorRampPalette(c("#2166AC", "white", "#CB181D"))(100),
  breaks        = seq(-abs_max, abs_max, length.out = 101),
  main          = paste0("ALS Upward Module — logFC ALS vs Ctrl\n(top ", n_heatmap, " genes by d120 logFC)"),
  filename      = "projects/Anderson-Suthar_collab/figures/als_up_module_heatmap.pdf",
  width = 5, height = 10
)

# ---- 7. Validate: overlap with WP ALS pathway gene set ----

als_pathway      <- gmtPathways("genesets/WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS.v2026.1.Hs.gmt")
als_pathway_genes <- unlist(als_pathway)

module_symbols <- als_up_module$symbol
overlap_genes  <- intersect(module_symbols, als_pathway_genes)

message("\n--- WP ALS Pathway Overlap ---")
message("Module size:          ", length(module_symbols))
message("Pathway size:         ", length(als_pathway_genes))
message("Overlap:              ", length(overlap_genes))
message("Overlap genes: ",        paste(overlap_genes, collapse = ", "))

# Fisher's exact test (one-sided: enrichment)
n_bg <- nrow(lfc_mat)
a <- length(overlap_genes)
b <- length(module_symbols) - a
c <- length(als_pathway_genes) - a
d <- n_bg - a - b - c

fisher_res <- fisher.test(matrix(c(a, b, c, d), 2, 2), alternative = "greater")
message("Fisher p = ", signif(fisher_res$p.value, 3),
        "  OR = ", signif(fisher_res$estimate, 3))

# ---- 8. GO enrichment of the module ----

entrez_map <- bitr(
  module_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

go_bp <- enrichGO(
  gene          = entrez_map$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

if (nrow(as.data.frame(go_bp)) > 0) {
  p_go <- dotplot(go_bp, showCategory = 20) +
    labs(title = "ALS Upward Module — GO Biological Process") +
    theme(axis.text.y = element_text(size = 9))

  ggsave(
    "projects/Anderson-Suthar_collab/figures/als_up_module_GO_BP.pdf",
    p_go, width = 8, height = 8
  )

  save_checkpoint(
    go_bp,
    name = "als_up_module_GO_BP",
    dir  = "projects/Anderson-Suthar_collab/data/r_objects"
  )
} else {
  message("No significant GO BP terms found at FDR < 0.05")
}
