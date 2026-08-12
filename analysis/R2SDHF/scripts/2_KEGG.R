# ---- Libraries ----
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(openxlsx)
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(tibble)
library(purrr)
# ---- Data Load ----
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD.rds")

#Ad Entrez Id column to data list

res_shrunk_list <- purrr::map(res_shrunk_list, function(df) {
  
  ensembl_keys <- rownames(df)
  ensembl_keys <- ensembl_keys[!is.na(ensembl_keys)]
  ensembl_keys <- as.character(ensembl_keys)
  
  map_df <- AnnotationDbi::select(
    x       = org.Hs.eg.db,
    keys    = unique(ensembl_keys),
    keytype = "ENSEMBL",
    columns = "ENTREZID"
  ) |>
    tibble::as_tibble() |>
    dplyr::rename(ensembl_id = ENSEMBL, entrez_id = ENTREZID) |>
    dplyr::filter(!is.na(ensembl_id)) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)
  
  df_tbl <- as.data.frame(df) |>
    tibble::rownames_to_column("ensembl_id") |>
    tibble::as_tibble()
  
  dplyr::left_join(df_tbl, map_df, by = "ensembl_id")
})



# ---- Kegg ----
#Extract signficant genes
p_col   <- "padj" 
p_cut   <- 0.05
lfc_cut <- 1

sigGenes_list <- purrr::imap(res_shrunk_list, function(df, contrast_name) {
  
  p_sym <- rlang::sym(p_col)
  
  df %>%
    tidyr::drop_na(entrez_id, log2FoldChange, !!p_sym) %>%
    dplyr::filter((!!p_sym) < p_cut, abs(log2FoldChange) > lfc_cut) %>%
    dplyr::pull(entrez_id) %>%
    unique()
})

keggRes_list <- purrr::imap(sigGenes_list, function(sigGenes, contrast_name) {
  
  if (length(sigGenes) == 0) return(NULL)
  
  clusterProfiler::enrichKEGG(
    gene     = sigGenes,
    organism = "hsa"
  )
})

# Kegg table
outdir  <- "projects/R2SDHF/data/kegg"
outfile <- file.path(outdir, "R2SDHF_kegg_results_noMAD.xlsx")

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


# ---- KEGG Pathview (ALS + Alzheimer) for each contrast ----
p_col   <- "padj"   # "pvalue" or "padj"
p_cut   <- 0.05
lfc_cut <- 1

pathways <- list(
  als       = "hsa05014",
  alzheimer = "hsa05010"
)

dir.create("projects/R2SDHF/figures/plasmo/pathways/als", recursive = TRUE, showWarnings = FALSE)
dir.create("projects/R2SDHF/figures/plasmo/pathways/alzheimer", recursive = TRUE, showWarnings = FALSE)

for (contrast_name in names(res_shrunk_list)) {
  
  df <- res_shrunk_list[[contrast_name]]
  
  logFC <- df$log2FoldChange
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
    
    outdir <- file.path("projects/R2SDHF/figures/plasmo/pathways", folder)
    
    withr::with_dir(outdir, {
      pathview::pathview(
        gene.data   = logFC,
        pathway.id  = pid,
        species     = "hsa",
        limit       = list(gene = 3, cpd = 1),
        out.suffix  = prefix,
        kegg.native = TRUE
      )
    })
  }
}


# ---- Top 5 pathways for each treatment ----

base_outdir <- "projects/R2SDHF/figures/plasmo/pathways/top5_FE"
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
  
  logFC <- df$log2FoldChange
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
        species     = "hsa",
        kegg.native = TRUE,
        out.suffix  = prefix,
        limit       = list(gene = 3, cpd = 1)
      )
    })
  }
}

# ---- GSEA ----
gmt_sym <- clusterProfiler::read.gmt(
  "genesets/WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS.v2026.1.Hs.gmt"
)

sym2ent <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gmt_sym$gene),
  keytype = "SYMBOL",
  columns = "ENTREZID"
)

gmt_ent <- gmt_sym |>
  dplyr::left_join(
    tibble::as_tibble(sym2ent) |>
      dplyr::rename(gene = SYMBOL, entrez_id = ENTREZID) |>
      dplyr::filter(!is.na(entrez_id)) |>
      dplyr::distinct(gene, .keep_all = TRUE),
    by = "gene"
  ) |>
  dplyr::filter(!is.na(entrez_id)) |>
  dplyr::transmute(term = term, gene = entrez_id)

ranked_list <- purrr::map(res_shrunk_list, function(df) {
  
  df2 <- df |>
    dplyr::filter(!is.na(entrez_id), !is.na(log2FoldChange)) |>
    dplyr::mutate(entrez_id = as.character(entrez_id)) |>
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
    dplyr::distinct(entrez_id, .keep_all = TRUE)
  
  geneList <- df2$log2FoldChange
  names(geneList) <- df2$entrez_id
  
  sort(geneList, decreasing = TRUE)
})

gsea_results <- purrr::map(ranked_list, function(geneList) {
  clusterProfiler::GSEA(
    geneList     = geneList,
    TERM2GENE    = gmt_ent,
    pvalueCutoff = 1
  )
})

gsea_tbl_list <- purrr::imap(gsea_results, function(res, contrast) {
  tibble::as_tibble(res@result) |>
    dplyr::mutate(contrast = contrast, .before = 1)
})


# ---- Table ----
# Assumes you already have:
# - res_shrunk_list (list of tibbles) with columns: entrez_id, log2FoldChange (or stat)
# - gmt_ent (TERM2GENE tibble/data.frame with columns: term, gene [ENTREZ IDs as character])

ranked_list <- purrr::map(res_shrunk_list, function(df) {
  
  df2 <- df |>
    dplyr::filter(!is.na(entrez_id), !is.na(log2FoldChange)) |>
    dplyr::mutate(entrez_id = as.character(entrez_id)) |>
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
    dplyr::distinct(entrez_id, .keep_all = TRUE)
  
  geneList <- df2$log2FoldChange
  names(geneList) <- df2$entrez_id
  
  sort(geneList, decreasing = TRUE)
})

gsea_results <- purrr::map(ranked_list, function(geneList) {
  clusterProfiler::GSEA(
    geneList = geneList,
    TERM2GENE = gmt_ent,
    pvalueCutoff = 1,
    minGSSize = 1
  )
})

gsea_tbl_all <- purrr::imap_dfr(gsea_results, function(res, contrast) {
  
  out <- tibble::as_tibble(res@result)
  
  if (nrow(out) == 0) {
    return(tibble::tibble(
      contrast = contrast,
      term = NA_character_,
      setSize = NA_integer_,
      NES = NA_real_,
      pvalue = NA_real_,
      p.adjust = NA_real_
    ))
  }
  
  out |>
    dplyr::transmute(
      contrast = contrast,
      term = ID,
      setSize,
      NES,
      pvalue,
      p.adjust
    )
})

gsea_tbl_all


write.table(
  gsea_tbl_all,
  file = "clipboard",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)











gsea_summary <- purrr::imap_dfr(gsea_results, function(res, contrast) {
  
  if (nrow(res@result) == 0) {
    return(
      tibble::tibble(
        contrast = contrast,
        NES = NA_real_,
        pvalue = NA_real_,
        p.adjust = NA_real_,
        size = NA_integer_,
        leading_edge = NA_character_
      )
    )
  }
  
  tibble::as_tibble(res@result) |>
    dplyr::transmute(
      contrast = contrast,
      NES,
      pvalue,
      p.adjust,
      size,
      leading_edge
    )
})

gsea_summary <- gsea_summary |>
  dplyr::mutate(
    NES = round(NES, 3),
    pvalue = signif(pvalue, 3),
    p.adjust = signif(p.adjust, 3)
  ) |>
  dplyr::arrange(p.adjust)
