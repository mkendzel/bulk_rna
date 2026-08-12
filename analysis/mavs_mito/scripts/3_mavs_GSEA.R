# =============================================
# 3 - GSEA (fgsea, MSigDB Hallmark)
# =============================================
# Loads the symbol-annotated limma-voom results from 2_mavs_LimmaVoom_DE.R
# and runs fgsea against the mouse Hallmark gene sets for:
#   - Spleen: WT vs MAVSKO at D7
#   - Brain:  WT vs MAVSKO at D7
#   - Brain vs Spleen at D7 (WT and MAVSKO, Spleen = reference)
# Ranking uses the moderated t statistic (see R/run_gsea.R).

# ---- Libraries ----
library(tidyverse)
library(fgsea)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

# ---- GSEA Spleen: WT vs MAVSKO D7 ----
res_voom_spleen_d0_d7_annotated <- load_checkpoint("res_voom_spleen_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects")

res_wt_vs_mavsko_d7 <- res_voom_spleen_d0_d7_annotated[["WT_vs_MAVSKO_D7"]]

fgsea_wt_vs_mavsko_d7 <- run_gsea(res_wt_vs_mavsko_d7, pathways)

# sanity checks
fgsea_wt_vs_mavsko_d7 %>% select(pathway, NES, pval, padj) %>% head(15)

fgsea_wt_vs_mavsko_d7 %>% filter(grepl("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|GLYCOLYSIS|HYPOXIA", pathway)) %>%
  select(pathway, NES, pval, padj)

sum(fgsea_wt_vs_mavsko_d7$padj < 0.05)

save_checkpoint(fgsea_wt_vs_mavsko_d7, "fgsea_wt_vs_mavsko_d7",
  notes = "GSEA results for WT vs MAVSKO at D7 in spleen. Pathways from MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(20:33)
)

# ---- GSEA Brain: WT vs MAVSKO D7 ----
res_voom_brain_d0_d7_annotated <- load_checkpoint("res_voom_brain_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects")

res_wt_vs_mavsko_d7_brain <- res_voom_brain_d0_d7_annotated[["WT_vs_MAVSKO_D7"]]

fgsea_wt_vs_mavsko_d7_brain <- run_gsea(res_wt_vs_mavsko_d7_brain, pathways)

# sanity checks
fgsea_wt_vs_mavsko_d7_brain %>% select(pathway, NES, pval, padj) %>% head(15)

fgsea_wt_vs_mavsko_d7_brain %>% filter(grepl("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|GLYCOLYSIS|HYPOXIA", pathway)) %>%
  select(pathway, NES, pval, padj)

sum(fgsea_wt_vs_mavsko_d7_brain$padj < 0.05)

save_checkpoint(fgsea_wt_vs_mavsko_d7_brain, "fgsea_wt_vs_mavsko_d7_brain",
  notes = "GSEA results for WT vs MAVSKO at D7 in brain. Pathways from MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(46:59)
)

# ---- GSEA: Brain vs Spleen D7 (Hallmark) ----
res_voom_d7_tissue_annotated <- load_checkpoint("res_voom_d7_tissue_annotated",
  dir = "projects/mavs_mito/data/r_objects")

fgsea_wt_brain_vs_spleen_d7     <- run_gsea(res_voom_d7_tissue_annotated[["WT_Brain_vs_Spleen"]], pathways)
fgsea_mavsko_brain_vs_spleen_d7 <- run_gsea(res_voom_d7_tissue_annotated[["MAVSKO_Brain_vs_Spleen"]], pathways)

# sanity checks (positive NES = enriched in Brain relative to Spleen)
fgsea_wt_brain_vs_spleen_d7 %>% select(pathway, NES, pval, padj) %>% head(15)
fgsea_mavsko_brain_vs_spleen_d7 %>% select(pathway, NES, pval, padj) %>% head(15)

sum(fgsea_wt_brain_vs_spleen_d7$padj < 0.05)
sum(fgsea_mavsko_brain_vs_spleen_d7$padj < 0.05)

save_checkpoint(fgsea_wt_brain_vs_spleen_d7, "fgsea_wt_brain_vs_spleen_d7",
  notes = "GSEA Brain vs Spleen at D7, WT (Spleen = ref). MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(72:82)
)

save_checkpoint(fgsea_mavsko_brain_vs_spleen_d7, "fgsea_mavsko_brain_vs_spleen_d7",
  notes = "GSEA Brain vs Spleen at D7, MAVSKO (Spleen = ref). MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(72:82)
)
