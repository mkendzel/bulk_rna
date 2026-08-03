# ---- Libraries ----
# Bioconductor first so dplyr masks S4Vectors::rename/select, not the reverse
library(org.Hs.eg.db)
library(edgeR)
library(limma)
library(dplyr)
library(stringr)
library(tibble)
library(purrr)
library(ggplot2)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("projects/WNV_ALS_R01_2026/scripts/0_config.R")

# ----- Import data ----
# Plasmidsaurus matrix: gene_id / gene_name / gene_biotype, then a <sample>_cpm and
# <sample>_count column per sample. Counts are fractional and get rounded below.
expr <- read.delim(
  expr_path,
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

# Kept for symbol fallback during annotation
vendor_genes <- expr[, intersect(c("gene_id", "gene_name", "gene_biotype"), colnames(expr))]

# ---- set up data frame ----
count_cols <- grep("_count$", colnames(expr), value = TRUE)

counts <- expr[, c("gene_id", count_cols)]
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

colnames(counts) <- sub("_count$", "", colnames(counts))

# ---- Rename plasmid codes to real sample names ----
# Match on the exact plasmid code, else on the trailing index (1-40), which
# survives a change in the vendor's column spelling.
rename_by_plasmid <- setNames(sample_map$sample, sample_map$plasmid)
rename_by_idx     <- setNames(sample_map$sample, as.character(sample_map$idx))

new_names <- vapply(colnames(counts), function(cn) {
  if (cn %in% names(rename_by_plasmid)) return(rename_by_plasmid[[cn]])
  idx <- sub("^.*_", "", cn)
  if (grepl("^[0-9]+$", idx) && idx %in% names(rename_by_idx)) return(rename_by_idx[[idx]])
  NA_character_
}, character(1), USE.NAMES = FALSE)

# List failures before stopping, so the fix is obvious
if (anyNA(new_names)) {
  message("Unmapped count columns:\n  ",
          paste(colnames(counts)[is.na(new_names)], collapse = "\n  "))
}
colnames(counts) <- new_names

stopifnot(
  !anyNA(colnames(counts)),
  !any(duplicated(colnames(counts))),
  length(setdiff(sample_map$sample, colnames(counts))) == 0
)

counts <- counts[, sample_map$sample, drop = FALSE]

# ---- QC For samples ----
qc_list <- lapply(colnames(counts), function(s) {
  x <- counts[, s, drop = TRUE]

  list(
    qc = data.frame(
      sample           = s,
      total_reads      = sum(x),
      detected_genes   = sum(x > 0),
      percent_of_total = sum(x) / sum(colSums(counts)) * 100
    ),
    gene_summary = summary(x)
  )
})
names(qc_list) <- colnames(counts)

qc_summary <- do.call(rbind, lapply(qc_list, `[[`, "qc"))
rownames(qc_summary) <- NULL

gene_summaries <- do.call(cbind, lapply(qc_list, function(z) z$gene_summary))
colnames(gene_summaries) <- names(qc_list)

qc_summary <- qc_summary |>
  dplyr::left_join(
    dplyr::select(sample_map, sample, line, stim, experiment),
    by = "sample"
  )

qc_summary

# Library size and gene detection per sample
ensure_dir(dir_fig, "qc")

qc_long <- qc_summary |>
  tidyr::pivot_longer(c(total_reads, detected_genes),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(
    sample = factor(sample, levels = sample_map$sample),
    metric = factor(metric, c("total_reads", "detected_genes"),
                    c("Total counts", "Detected genes"))
  )

qc_plot <- ggplot(qc_long, aes(x = sample, y = value, fill = line)) +
  geom_col() +
  facet_grid(metric ~ experiment, scales = "free", space = "free_x") +
  scale_fill_manual(values = line_cols) +
  scale_y_continuous(labels = scales::comma) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(dir_fig, "qc", "library_size_and_detection.png"),
       qc_plot, width = 12, height = 6, dpi = 300)

# EDIT: drop low-count samples listed in qc_summary above
drop_samples <- character(0)
counts <- counts[, setdiff(colnames(counts), drop_samples), drop = FALSE]

# Recompute replication from the samples that survived, since scripts 2-3 read
# min_n to flag unreplicated contrasts
cond_n_kept <- sample_map |>
  dplyr::filter(sample %in% colnames(counts)) |>
  dplyr::count(condition, .drop = FALSE, name = "n") |>
  tibble::deframe()

contrast_registry$min_n <- vapply(
  contrast_registry$expr,
  function(e) {
    g <- intersect(unlist(strsplit(e, "[^A-Za-z0-9_]+")), names(cond_n_kept))
    min(cond_n_kept[g])
  },
  integer(1)
)

if (length(drop_samples) > 0) {
  message("Dropped ", length(drop_samples), " sample(s). Groups now below 2 replicates: ",
          paste(names(cond_n_kept)[cond_n_kept < 2], collapse = ", "))
}

# Save/load checkpoint
save_checkpoint(counts, "counts_filtered", dir = dir_rds)
# counts <- load_checkpoint("counts_filtered", dir = dir_rds)

save_checkpoint(contrast_registry, "contrast_registry", dir = dir_rds)

# ---- limma-voom fit, one experiment at a time ----
# Experiments are fit separately, so library sizes and dispersion are never shared
# between the spinal and cortical batches. The cell-means design (~ 0 + condition)
# makeContrasts express the interaction terms.
fit_experiment <- function(exp_name, counts, registry) {

  sm <- sample_map |>
    dplyr::filter(experiment == exp_name, sample %in% colnames(counts)) |>
    dplyr::mutate(condition = droplevels(condition))

  mat <- round(as.matrix(counts[, sm$sample, drop = FALSE]))

  design <- model.matrix(~ 0 + condition, data = sm)
  colnames(design) <- levels(sm$condition)

  dge <- edgeR::DGEList(counts = mat, group = sm$condition)
  dge$samples <- cbind(
    dge$samples,
    as.data.frame(sm[, c("sample", "line", "stim", "rep", "experiment")])
  )

  keep <- edgeR::filterByExpr(dge, design = design)
  message(exp_name, ": keeping ", sum(keep), " of ", length(keep), " genes")
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- edgeR::calcNormFactors(dge)

  png(file.path(dir_fig, "qc", paste0("voom_mean_variance_", exp_name, ".png")),
      width = 7, height = 5, units = "in", res = 300)
  v <- limma::voom(dge, design, plot = TRUE)
  dev.off()

  fit <- limma::lmFit(v, design)

  reg <- registry[registry$experiment == exp_name, ]

  # Drop contrasts whose groups lost all samples during QC
  has_groups <- vapply(reg$expr, function(e) {
    all(intersect(unlist(strsplit(e, "[^A-Za-z0-9_]+")), levels(sample_map$condition))
        %in% colnames(design))
  }, logical(1))

  if (any(!has_groups)) {
    message(exp_name, ": skipping ", sum(!has_groups), " contrast(s) with no remaining samples: ",
            paste(reg$name[!has_groups], collapse = ", "))
    reg <- reg[has_groups, ]
  }

  cm <- limma::makeContrasts(contrasts = reg$expr, levels = design)
  colnames(cm) <- reg$name

  fit2 <- limma::eBayes(limma::contrasts.fit(fit, cm))

  tt <- setNames(
    lapply(reg$name, function(cn) {
      limma::topTable(fit2, coef = cn, number = Inf, sort.by = "none")
    }),
    reg$name
  )

  list(sample_meta = sm, design = design, dge = dge,
       voom = v, fit = fit2, contrast_matrix = cm, toptables = tt)
}

fits <- setNames(
  lapply(EXPERIMENTS, fit_experiment, counts = counts, registry = contrast_registry),
  EXPERIMENTS
)

# Expect 31 x 12 (spinal) and 9 x 3 (cort)
lapply(fits, function(f) dim(f$design))
lapply(fits, function(f) length(f$toptables))

# ---- Map gene IDs ----
all_ensembl <- fits |>
  purrr::map(~ purrr::map(.x$toptables, rownames)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Strip version suffix (ENSG00000123456.7 -> ENSG00000123456)
all_ensembl_clean <- unique(sub("\\..*$", "", all_ensembl))

# Human: org.Hs.eg.db. Mouse: org.Mm.eg.db
gene_map <- suppressMessages(AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensembl_clean,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)) |>
  tibble::as_tibble() |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID
  )

# Fall back to the vendor gene name where org.Hs.eg.db has no symbol, so fewer
# genes are unrankable in script 2's GSEA
if ("gene_name" %in% colnames(vendor_genes)) {
  vendor_sym <- vendor_genes |>
    dplyr::transmute(
      ensembl_id = sub("\\..*$", "", gene_id),
      gene_name  = dplyr::na_if(gene_name, "")
    ) |>
    dplyr::distinct(ensembl_id, .keep_all = TRUE)

  gene_map <- gene_map |>
    dplyr::left_join(vendor_sym, by = "ensembl_id") |>
    dplyr::mutate(SYMBOL = dplyr::coalesce(SYMBOL, gene_name)) |>
    dplyr::select(-gene_name)
}

# ---- Append identifiers to every contrast ----
annotate_toptables <- function(tt_list, gene_map) {
  purrr::map(tt_list, function(df) {
    df |>
      as.data.frame() |>
      tibble::rownames_to_column("ensembl_id") |>
      dplyr::mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
      dplyr::left_join(gene_map, by = "ensembl_id") |>
      tibble::as_tibble()
  })
}

tt_annotated <- purrr::map(fits, ~ annotate_toptables(.x$toptables, gene_map))

# ---- Save checkpoints ----
# Downstream columns are limma's: logFC, AveExpr, t, P.Value, adj.P.Val, B.
# Expression values are voom log2-CPM (v$E).
for (e in EXPERIMENTS) {
  save_checkpoint(fits[[e]]$dge,         paste0("dge_", e),              dir = dir_rds)
  save_checkpoint(fits[[e]]$voom,        paste0("voom_", e),             dir = dir_rds)
  save_checkpoint(fits[[e]]$fit,         paste0("fit_", e),              dir = dir_rds)
  save_checkpoint(fits[[e]]$voom$E,      paste0("logcpm_", e),           dir = dir_rds)
  save_checkpoint(fits[[e]]$sample_meta, paste0("sample_meta_", e),      dir = dir_rds)
  save_checkpoint(tt_annotated[[e]],     paste0("tt_", e, "_annotated"), dir = dir_rds)
}

# tt_annotated <- setNames(
#   lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)),
#   EXPERIMENTS
# )

# ---- DE counts per contrast ----
de_counts <- purrr::imap_dfr(tt_annotated, function(tts, exp_name) {
  purrr::imap_dfr(tts, function(df, cn) {
    tibble::tibble(
      contrast = cn,
      up       = sum(df$adj.P.Val < P_CUT &  df$logFC >  LFC_CUT, na.rm = TRUE),
      down     = sum(df$adj.P.Val < P_CUT &  df$logFC < -LFC_CUT, na.rm = TRUE)
    )
  })
}) |>
  dplyr::left_join(contrast_registry, by = c("contrast" = "name")) |>
  dplyr::select(experiment, type, contrast, label, min_n, up, down) |>
  dplyr::arrange(experiment, type, contrast)

de_counts

ensure_dir(dir_res)
write.table(de_counts, file.path(dir_res, "de_gene_counts.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
