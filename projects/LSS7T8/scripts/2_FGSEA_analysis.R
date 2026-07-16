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