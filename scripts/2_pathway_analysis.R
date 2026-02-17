# ---- Libraries ----
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)

# install pathview if not already installed from Bioconductor
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}


# ---- Add Identifier Mapping ----
res_shrunk_list <- purrr::map(res_shrunk_list, function(df) {
  
  map_df <- AnnotationDbi::select(
    x       = org.Hs.eg.db,
    keys    = unique(df$ensembl_id),
    keytype = "ENSEMBL",
    columns = c("ENTREZID")
  ) %>%
    as_tibble() %>%
    dplyr::rename(ensembl_id = 1) %>%
    dplyr::filter(!is.na(ensembl_id)) %>%
    dplyr::distinct(ensembl_id, .keep_all = TRUE)
  
  if ("ENTREZID" %in% names(map_df)) {
    map_df <- map_df %>% dplyr::rename(entrez_id = ENTREZID)
  } else {
    map_df <- map_df %>% dplyr::rename(entrez_id = 2)
  }
  
  df %>%
    dplyr::left_join(map_df, by = "ensembl_id")
})


# ---- Create df -----
p_col   <- "padj"   # "pvalue" or "padj"
p_cut   <- 0.05
lfc_cut <- 1.5

pathway_ids <- c("hsa05020")

dir.create("figures/plasmo/hsa05020/pathways", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/plasmo/pathways/hsa05020/kegg_tables", recursive = TRUE, showWarnings = FALSE)

kegg_tbl_list <- purrr::imap(res_shrunk_list, function(df, contrast_name) {
  
  p_sym <- rlang::sym(p_col)
  
  sigGenes <- df %>%
    tidyr::drop_na(entrez_id, log2FoldChange, !!p_sym) %>%
    dplyr::filter((!!p_sym) < p_cut, abs(log2FoldChange) > lfc_cut) %>%
    dplyr::pull(entrez_id) %>%
    unique()
  
  keggRes <- clusterProfiler::enrichKEGG(gene = sigGenes, organism = "hsa")
  
  kegg_tbl <- tibble::as_tibble(keggRes) %>%
    dplyr::mutate(
      contrast = contrast_name,
      p_col = p_col,
      p_cutoff = p_cut,
      lfc_cutoff = lfc_cut,
      .before = 1
    )
  
  readr::write_csv(
    kegg_tbl,
    file = file.path("figures/plasmo/pathways/hsa05020/kegg_tables", paste0(contrast_name, "_kegg.csv"))
  )
  
  logFC <- df$log2FoldChange
  names(logFC) <- df$entrez_id
  logFC <- logFC[!is.na(names(logFC))]
  
  out_dir <- file.path("figures/plasmo/hsa05020/pathways", contrast_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  withr::with_dir(out_dir, {
    purrr::walk(pathway_ids, ~{
      pathview::pathview(
        gene.data  = logFC,
        pathway.id = .x,
        species    = "hsa",
        limit      = list(gene = 2, cpd = 1),
        out.suffix = contrast_name,
        kegg.dir   = "."
      )
    })
  })
  
  kegg_tbl
})

