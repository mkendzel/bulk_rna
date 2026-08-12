# ---- Libraries ----
library(tidyverse)
library(clusterProfiler)
library(pathview)
library(openxlsx)
library(org.Mm.eg.db)
# ---- Load helper functions/Data Load ----
invisible(sapply(list.files("R", full.names = TRUE), source))

res_shrunk_list <- load_checkpoint("res_voom_spleen_d0_d7_annotated", dir = "projects/mavs_mito/data/r_objects")

# ---- Add mouse Entrez IDs (KEGG needs Entrez, not Ensembl) ----
res_shrunk_list <- lapply(res_shrunk_list, function(df) {
  df$entrez_id <- mapIds(
    org.Mm.eg.db,
    keys      = df$ensembl_id,
    keytype   = "ENSEMBL",
    column    = "ENTREZID",
    multiVals = "first"
  )
  df
})

# ---- Kegg ----
#Extract signficant genes based on lfc
p_col   <- "adj.P.Val"
p_cut   <- 0.05
lfc_cut <- 1

sigGenes_list <- purrr::imap(res_shrunk_list, function(df, contrast_name) {
  
  p_sym <- rlang::sym(p_col)
  
  df %>%
    tidyr::drop_na(entrez_id, logFC, !!p_sym) %>%
    dplyr::filter((!!p_sym) < p_cut, abs(logFC) > lfc_cut) %>%
    dplyr::pull(entrez_id) %>%
    unique()
})

# Run KEGG enrichment for each contrast
keggRes_list <- purrr::imap(sigGenes_list, function(sigGenes, contrast_name) {

  if (length(sigGenes) == 0) return(NULL)

  clusterProfiler::enrichKEGG(
    gene     = sigGenes,
    organism = "mmu"  # Mouse: "mmu"
  )
})

# Kegg table save
outdir  <- "projects/mavs_mito/data/kegg"
outfile <- file.path(outdir, "kegg_results.xlsx")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

wb <- openxlsx::createWorkbook()

purrr::imap(keggRes_list, function(res, contrast_name) {
  
  if (is.null(res)) return(NULL)
  
  tbl <- tibble::as_tibble(res)
  
  sheet_name <- gsub("[^A-Za-z0-9]", "_", contrast_name)
  sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name limit
  
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, tbl)
})

openxlsx::saveWorkbook(wb, outfile, overwrite = TRUE)


# ---- KEGG Pathview for selected pathways, for each contrast ----
p_col   <- "adj.P.Val"   # limma-voom column
p_cut   <- 0.05
lfc_cut <- 1

# EDIT: pick pathway IDs (use clusterProfiler::search_kegg_organism()).
# Names below become the output subfolder names.
pathways <- list(
  T_cell_receptor_signaling = "mmu04110"
)

for (folder in names(pathways)) {
  dir.create(file.path("projects/mavs_mito/figures/pathways", folder), recursive = TRUE, showWarnings = FALSE)
}

for (contrast_name in names(res_shrunk_list)) {

  df <- res_shrunk_list[[contrast_name]]

  logFC <- df$logFC
  names(logFC) <- df$entrez_id
  logFC <- logFC[!is.na(names(logFC))]
  logFC <- logFC[!is.na(logFC)]

  for (folder in names(pathways)) {

    pid <- pathways[[folder]]

    safe_contrast <- gsub("[^A-Za-z0-9_]+", "_", contrast_name)
    prefix <- paste0(
      safe_contrast,
      "__", p_col, "lt", p_cut,
      "__lfcgt", lfc_cut,
      "__", pid, "_"
    )

    outdir <- file.path("projects/mavs_mito/figures/pathways", folder)

    withr::with_dir(outdir, {
      pathview::pathview(
        gene.data   = logFC,
        pathway.id  = pid,
        species     = "mmu",  # Mouse
        limit       = list(gene = 3, cpd = 1),
        out.suffix  = prefix,
        kegg.native = TRUE
      )
    })
  }
}


# ---- Top 5 pathways for each treatment ----

base_outdir <- "projects/mavs_mito/figures/pathways/top5_FE"
dir.create(base_outdir, recursive = TRUE, showWarnings = FALSE)

ratio_to_num <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(p) as.numeric(p[1]) / as.numeric(p[2]), numeric(1))
}

for (contrast_name in names(keggRes_list)) {
  
  keggRes <- keggRes_list[[contrast_name]]
  if (is.null(keggRes)) next
  
  kegg_tbl <- tibble::as_tibble(keggRes)
  
  if (!all(c("ID", "GeneRatio", "BgRatio") %in% names(kegg_tbl))) next
  
  top5_ids <- kegg_tbl |>
    dplyr::mutate(
      FoldEnrichment = ratio_to_num(GeneRatio) / ratio_to_num(BgRatio)
    ) |>
    dplyr::arrange(dplyr::desc(FoldEnrichment)) |>
    dplyr::slice_head(n = 5) |>
    dplyr::pull(ID)
  
  df <- res_shrunk_list[[contrast_name]]

  logFC <- df$logFC
  names(logFC) <- df$entrez_id
  logFC <- logFC[!is.na(names(logFC))]
  logFC <- logFC[!is.na(logFC)]

  contrast_dir <- file.path(base_outdir, gsub("[^A-Za-z0-9_]+", "_", contrast_name))
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (pid in top5_ids) {

    prefix <- paste0(
      gsub("[^A-Za-z0-9_]+", "_", contrast_name),
      "__top5FE__",
      pid,
      "_"
    )

    withr::with_dir(contrast_dir, {
      pathview::pathview(
        gene.data   = logFC,
        pathway.id  = pid,
        species     = "mmu",  # Mouse
        kegg.native = TRUE,
        out.suffix  = prefix,
        limit       = list(gene = 3, cpd = 1)
      )
    })
  }
}
