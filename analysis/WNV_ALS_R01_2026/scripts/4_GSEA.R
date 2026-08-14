# GSEA on the limma contrasts from 3_DE_limma.R. Compute only - the Hallmark
# bubble and ALS bar figures are drawn in 5_Graphs.R from the checkpoints below.
# Run from the repo root.

# ---- Libraries ----
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Dataset ----
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

tt <- setNames(
  lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)),
  EXPERIMENTS
)

contrast_registry <- load_checkpoint("contrast_registry", dir = dir_rds)

ensure_dir(dir_res)

# ---- Hallmark GSEA ----
# run_gsea() (R/run_gsea.R) ranks by moderated t, de-duplicating symbols on
# highest AveExpr
gsea_results <- purrr::map(tt, function(tts) {
  purrr::imap(tts, function(df, contrast_name) {
    message("fgsea: ", contrast_name)
    run_gsea(as.data.frame(df), hallmark_sets, gene_col = "SYMBOL")
  })
})

# leadingEdge is a list column; collapse it for writing to disk
flatten_le <- function(x) {
  if (!"leadingEdge" %in% names(x)) return(x)
  x$leadingEdge <- vapply(x$leadingEdge, paste, character(1), collapse = ";")
  x
}

gsea_tbl_all <- purrr::imap_dfr(gsea_results, function(gr, exp_name) {
  purrr::imap_dfr(gr, function(res, contrast_name) {
    tibble::as_tibble(flatten_le(as.data.frame(res))) |>
      dplyr::mutate(experiment = exp_name, contrast = contrast_name, .before = 1)
  })
}) |>
  dplyr::left_join(
    dplyr::select(contrast_registry, name, type, line, ref_line, stim, label, min_n),
    by = c("contrast" = "name")
  ) |>
  dplyr::arrange(experiment, type, contrast, padj)

write.table(gsea_tbl_all, file.path(dir_res, "gsea_hallmark_all.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# 5_Graphs.R draws the bubble plot from this; the RDS keeps `type` a factor
save_checkpoint(gsea_tbl_all, "gsea_hallmark", dir = dir_rds)

# ---- Custom geneset GSEA (clusterProfiler) ----
# EDIT: genesets/HP_ALZHEIMER_DISEASE.v2026.1.Hs.gmt is the other set available
custom_gmt <- "genesets/WP_AMYOTROPHIC_LATERAL_SCLEROSIS_ALS.v2026.1.Hs.gmt"

gmt_sym <- clusterProfiler::read.gmt(custom_gmt)

sym2ent <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = unique(gmt_sym$gene),
  keytype = "SYMBOL",
  columns = "ENTREZID"
))

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

# Rank on moderated t, matching the Hallmark ranking above
ranked_list_entrez <- purrr::map(tt, function(tts) {
  purrr::map(tts, function(df) {

    df2 <- df |>
      dplyr::filter(!is.na(entrez_id), !is.na(t)) |>
      dplyr::mutate(entrez_id = as.character(entrez_id)) |>
      dplyr::arrange(dplyr::desc(abs(t))) |>
      dplyr::distinct(entrez_id, .keep_all = TRUE)

    geneList <- df2$t
    names(geneList) <- df2$entrez_id

    sort(geneList, decreasing = TRUE)
  })
})

# pvalueCutoff/minGSSize = 1 keeps every result for downstream filtering
gsea_results_custom <- purrr::map(ranked_list_entrez, function(rl) {
  purrr::imap(rl, function(geneList, contrast_name) {
    message("custom GSEA: ", contrast_name)
    suppressWarnings(clusterProfiler::GSEA(
      geneList     = geneList,
      TERM2GENE    = gmt_ent,
      pvalueCutoff = 1,
      minGSSize    = 1
    ))
  })
})

gsea_tbl_all_custom <- purrr::imap_dfr(gsea_results_custom, function(gr, exp_name) {
  purrr::imap_dfr(gr, function(res, contrast_name) {

    out <- tibble::as_tibble(res@result)

    if (nrow(out) == 0) {
      return(tibble::tibble(
        experiment = exp_name,
        contrast   = contrast_name,
        term       = NA_character_,
        setSize    = NA_integer_,
        NES        = NA_real_,
        pvalue     = NA_real_,
        p.adjust   = NA_real_
      ))
    }

    out |>
      dplyr::transmute(
        experiment = exp_name,
        contrast   = contrast_name,
        term       = ID,
        setSize,
        NES,
        pvalue,
        p.adjust
      )
  })
}) |>
  dplyr::left_join(
    dplyr::select(contrast_registry, name, type, label, min_n),
    by = c("contrast" = "name")
  ) |>
  dplyr::mutate(
    NES      = round(NES, 3),
    pvalue   = signif(pvalue, 3),
    p.adjust = signif(p.adjust, 3)
  ) |>
  dplyr::arrange(experiment, p.adjust)

gsea_tbl_all_custom

write.table(
  gsea_tbl_all_custom,
  file = file.path(dir_res, "gsea_als_all.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# 5_Graphs.R draws the NES bar chart from this
save_checkpoint(gsea_tbl_all_custom, "gsea_als", dir = dir_rds)
