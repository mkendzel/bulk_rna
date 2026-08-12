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
# ---- KEGG Toll-like receptor signaling: dopc vs cl ----
# =================================================
# The MSigDB KEGG TLR set ships as human symbols; these data are mouse.
# Title-case the human symbols to the mouse convention (Akt1, Myd88, ...).

library(KEGGREST)
library(org.Mm.eg.db)
library(AnnotationDbi)

# Fetch KEGG gene entries linked to the mouse TLR pathway (values like "mmu:21898")
tlr_kegg_ids <- keggLink("mmu", "path:mmu04620")
tlr_entrez   <- sub("^mmu:", "", unname(tlr_kegg_ids))

# Map Entrez IDs to mouse gene symbols
tlr_symbols <- mapIds(
  org.Mm.eg.db,
  keys    = unique(tlr_entrez),
  column  = "SYMBOL",
  keytype = "ENTREZID"
)
tlr_symbols <- unique(na.omit(unname(tlr_symbols)))

# Assemble a named list matching fgsea::gmtPathways() output
kegg_tlr <- list(KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY = tlr_symbols)

fgsea_tlr_dopc_vs_cl_liver  <- run_gsea(res_dopc_vs_cl_liver,  kegg_tlr)
fgsea_tlr_dopc_vs_cl_spleen <- run_gsea(res_dopc_vs_cl_spleen, kegg_tlr)

save_checkpoint(fgsea_tlr_dopc_vs_cl_liver, "fgsea_tlr_dopc_vs_cl_liver",
                dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA of KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY (mouse symbols) for dopc vs cl in liver, voom limma results")

save_checkpoint(fgsea_tlr_dopc_vs_cl_spleen, "fgsea_tlr_dopc_vs_cl_spleen",
                dir = "projects/LSS7T8/data/r_objects",
                notes = "FgSEA of KEGG_TOLL_LIKE_RECEPTOR_SIGNALING_PATHWAY (mouse symbols) for dopc vs cl in spleen, voom limma results")

# ---- Print KEGG TLR results ----
print_tlr <- function(fg, tt, tissue) {
  cat("\n==== KEGG TLR signaling | dopc vs cl |", tissue, "====\n")
  cat("Set genes:", length(kegg_tlr[[1]]),
      "| detected in results:", sum(unique(tt$SYMBOL) %in% kegg_tlr[[1]]), "\n")
  print(as.data.frame(fg[, c("pathway", "pval", "padj", "NES", "size")]))
  cat("Direction:", ifelse(fg$NES[1] > 0, "up in dopc", "up in cl"), "\n")
  cat("Leading edge (", length(fg$leadingEdge[[1]]), "):\n  ",
      paste(fg$leadingEdge[[1]], collapse = ", "), "\n", sep = "")
}

print_tlr(fgsea_tlr_dopc_vs_cl_liver,  res_dopc_vs_cl_liver,  "liver")
print_tlr(fgsea_tlr_dopc_vs_cl_spleen, res_dopc_vs_cl_spleen, "spleen")