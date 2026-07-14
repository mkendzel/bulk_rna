# ---- Libraries ----
library(fgsea)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyr)
library(dplyr)
library(openxlsx)
library(purrr)
library(tibble)
library(ggplot2)
# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

# ---- Dataset ----
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)

res_shrunk_list <- load_checkpoint("res_shrunk_list_annotated", dir = "projects/PROJECT_NAME/data/r_objects")

# ---- order ----
ranked_list <- res_shrunk_list %>%
  purrr::map(~{
    
    df <- as.data.frame(.x)
    df <- df[!is.na(df$SYMBOL), ]
    
    ranks_logfc <- df$log2FoldChange
    names(ranks_logfc) <- df$SYMBOL
    
    ranks_logfc <- ranks_logfc[is.finite(ranks_logfc)]
    ranks_logfc <- ranks_logfc[!duplicated(names(ranks_logfc))]
    ranks_logfc <- sort(ranks_logfc, decreasing = TRUE)
    
    ranks_logfc
  })


gsea_results <- map(
  ranked_list,
  ~ fgsea(
    pathways = hallmark_sets,
    stats    = .x
  )
)

# ---- GSEA Results ----

gsea_tbl_list <- imap(gsea_results, ~{
  as_tibble(.x) %>%
    arrange(padj) %>%
    mutate(contrast = .y, .before = 1)
})

# ---- Bubble graph ----
# Parses names like condition_TreatmentA_1_vs_Control_0
# EDIT: assumes <treatment>_<timepoint>_vs_<control> naming from script 1.
# Adjust if your design has no time axis.
parse_contrast <- function(comparison) {
  x <- sub("^condition_", "", comparison)
  parts <- strsplit(x, "_vs_", fixed = TRUE)[[1]]
  treat <- parts[1]
  ctrl  <- if (length(parts) > 1) parts[2] else NA_character_
  
  tr <- strsplit(treat, "_", fixed = TRUE)[[1]]
  treatment <- tr[1]
  time      <- suppressWarnings(as.integer(tr[2]))
  
  list(
    treat = treat,
    ctrl = ctrl,
    treatment = treatment,
    time = time
  )
}

# Combine all contrasts into one table
df_plot <- purrr::imap_dfr(gsea_results, ~{
  as_tibble(.x) %>%
    mutate(comparison = .y, .before = 1)
})

# Parse labels for x axis
ps <- lapply(df_plot$comparison, parse_contrast)
df_plot$treat     <- vapply(ps, `[[`, character(1), "treat")
df_plot$ctrl       <- vapply(ps, `[[`, character(1), "ctrl")
df_plot$treatment  <- vapply(ps, `[[`, character(1), "treatment")
df_plot$time       <- vapply(ps, `[[`, integer(1),   "time")

# Clean pathway labels (Hallmark)
df_plot$pathway <- sub("^HALLMARK_", "", df_plot$pathway)

# Keep significant pathways only
# EDIT: match the significance cutoff used elsewhere
df_plot <- df_plot %>%
  filter(!is.na(padj), padj <= 0.05, !is.na(NES), !is.na(pathway))

# Order x-axis by time, then treatment (within time)
x_levels <- df_plot %>%
  distinct(treat, treatment, time) %>%
  arrange(time, treatment) %>%
  pull(treat)

df_plot$treat <- factor(df_plot$treat, levels = x_levels)

# Order y-axis by mean NES
df_plot$pathway <- with(df_plot, reorder(pathway, NES, mean))

# Bubble plot
GSEA_bubble <- ggplot(
  df_plot,
  aes(x = treat, y = pathway, size = -log10(padj), fill = NES)
) +
  geom_point(shape = 21, color = "black") +
  scale_size(name = "-log10(padj)", range = c(4, 10)) +
  ylab("Gene Set") +
  scale_x_discrete(position = "top") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text.y = element_text(colour = "black"),
    axis.text.x = element_text(colour = "black", angle = 45, hjust = 0),
    axis.title.x = element_blank()
  ) +
  scale_fill_distiller(palette = "Spectral")

GSEA_bubble

# ---- Custom geneset GSEA (clusterProfiler) ----
# EDIT: pick a geneset relevant to your project's hypothesis/biology
gmt_sym <- clusterProfiler::read.gmt(
  "genesets/<GENESET_FILE>.gmt"
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

ranked_list_entrez <- purrr::map(res_shrunk_list, function(df) {

  df2 <- df |>
    dplyr::filter(!is.na(entrez_id), !is.na(log2FoldChange)) |>
    dplyr::mutate(entrez_id = as.character(entrez_id)) |>
    dplyr::arrange(dplyr::desc(abs(log2FoldChange))) |>
    dplyr::distinct(entrez_id, .keep_all = TRUE)

  geneList <- df2$log2FoldChange
  names(geneList) <- df2$entrez_id

  sort(geneList, decreasing = TRUE)
})

# pvalueCutoff/minGSSize = 1 keeps every result for downstream filtering
gsea_results_custom <- purrr::map(ranked_list_entrez, function(geneList) {
  clusterProfiler::GSEA(
    geneList     = geneList,
    TERM2GENE    = gmt_ent,
    pvalueCutoff = 1,
    minGSSize    = 1
  )
})

gsea_tbl_list_custom <- purrr::imap(gsea_results_custom, function(res, contrast) {
  tibble::as_tibble(res@result) |>
    dplyr::mutate(contrast = contrast, .before = 1)
})

# ---- Custom geneset GSEA: table ----
gsea_tbl_all_custom <- purrr::imap_dfr(gsea_results_custom, function(res, contrast) {

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

gsea_tbl_all_custom

dir.create("projects/PROJECT_NAME/results", recursive = TRUE, showWarnings = FALSE)
write.table(
  gsea_tbl_all_custom,
  file = "projects/PROJECT_NAME/results/gsea_summary.tsv",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

gsea_summary_custom <- purrr::imap_dfr(gsea_results_custom, function(res, contrast) {

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

gsea_summary_custom <- gsea_summary_custom |>
  dplyr::mutate(
    NES = round(NES, 3),
    pvalue = signif(pvalue, 3),
    p.adjust = signif(p.adjust, 3)
  ) |>
  dplyr::arrange(p.adjust)

# ---- Other organoids ----
