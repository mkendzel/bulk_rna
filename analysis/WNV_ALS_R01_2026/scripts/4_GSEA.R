# ---- Libraries ----
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(tidyr)
library(dplyr)
library(purrr)
library(tibble)
library(ggplot2)

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
ensure_dir(dir_fig, "gsea")

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

# ---- Bubble graph ----
make_bubble <- function(exp_name) {

  df_plot <- gsea_tbl_all |>
    dplyr::filter(experiment == exp_name,
                  !is.na(padj), padj <= P_CUT,
                  !is.na(NES), !is.na(pathway))

  if (nrow(df_plot) == 0) {
    message("No significant Hallmark pathways for ", exp_name)
    return(NULL)
  }

  df_plot$pathway <- sub("^HALLMARK_", "", df_plot$pathway)

  # x-axis follows registry order: within_line -> genotype -> interaction
  x_levels <- contrast_registry |>
    dplyr::filter(experiment == exp_name) |>
    dplyr::arrange(type, line, ref_line, stim) |>
    dplyr::pull(label)

  df_plot$label   <- factor(df_plot$label, levels = unique(x_levels))
  df_plot$pathway <- with(df_plot, reorder(pathway, NES, mean))
  # Unreplicated groups drawn as triangles
  df_plot$replication <- ifelse(df_plot$min_n == 1, "n = 1", "replicated")

  ggplot(df_plot, aes(x = label, y = pathway,
                      size = -log10(padj), fill = NES, shape = replication)) +
    geom_point(colour = "black") +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    scale_size(name = "-log10(padj)", range = c(3, 9)) +
    scale_shape_manual(values = c("n = 1" = 24, "replicated" = 21), drop = FALSE) +
    ylab("Gene Set") +
    scale_x_discrete(position = "top") +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(colour = "black"),
      axis.text.x = element_text(colour = "black", angle = 45, hjust = 0),
      axis.title.x = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = NA)
    ) +
    scale_fill_distiller(palette = "Spectral")
}

for (e in EXPERIMENTS) {
  p <- make_bubble(e)
  if (is.null(p)) next
  n_x <- dplyr::n_distinct(dplyr::filter(contrast_registry, experiment == e)$label)
  ggsave(file.path(dir_fig, "gsea", paste0("hallmark_bubble_", e, ".png")),
         p, width = max(7, 2.2 + n_x * 0.9), height = 9, dpi = 300, limitsize = FALSE)
}

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

# ---- ALS geneset NES per contrast ----
sig_label <- paste0("padj < ", P_CUT)

als_plot_df <- gsea_tbl_all_custom |>
  dplyr::filter(!is.na(NES))

if (nrow(als_plot_df) > 0) {

  als_plot_df <- als_plot_df |>
    dplyr::mutate(
      label = factor(label, levels = contrast_registry$label[
        order(contrast_registry$type, contrast_registry$line,
              contrast_registry$ref_line, contrast_registry$stim)
      ]),
      significant = ifelse(p.adjust < P_CUT, sig_label, "ns")
    )

  als_bar <- ggplot(als_plot_df, aes(x = label, y = NES, fill = significant)) +
    geom_col(colour = "black", width = 0.7) +
    facet_grid(~ type, scales = "free_x", space = "free_x") +
    geom_hline(yintercept = 0) +
    scale_fill_manual(values = setNames(c("#D62728", "grey80"), c(sig_label, "ns"))) +
    labs(x = NULL, y = "NES", title = "WP: Amyotrophic lateral sclerosis") +
    theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, colour = "black"),
      strip.background = element_rect(fill = "grey92", colour = NA)
    )

  ggsave(file.path(dir_fig, "gsea", "als_geneset_NES.png"),
         als_bar, width = 11, height = 5, dpi = 300)
}
