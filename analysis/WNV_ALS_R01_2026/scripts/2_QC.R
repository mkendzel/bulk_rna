# QC on the counts checkpoint from 1_df_count.R. Writes results/qc_report.md and
# figures/qc/. Nothing is filtered and no checkpoint is written.
# Run from the repo root.

# ---- Libraries ----
library(DESeq2)
library(RNAseqQC)
library(dplyr)
library(tibble)
library(ggplot2)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Load from 1_df_count.R ----
counts       <- load_checkpoint("annotated_counts", dir = dir_rds)
vendor_genes <- load_checkpoint("vendor_genes", dir = dir_rds)

# ---- QC thresholds ----
# Every check reads warn/fail from this table. Delete a row to drop that check.
#   direction  "higher" = larger is better, "lower" = smaller is better
#   warn/fail  absolute values in the metric's own unit
#   unit       "count" | "pct" | "ratio", formatting only
QC_THRESHOLDS <- tibble::tribble(
  ~metric,                ~label,                    ~category,     ~direction, ~warn,   ~fail,   ~unit,
  "uniquely_mapped",      "Uniquely mapped",         "alignment",   "higher",   10e6,    5e6,     "count",
  "dedup_unique_reads",   "Unique after dedup",      "alignment",   "higher",   6e6,     4e6,     "count",
  "mapping_rate",         "Mapping rate",            "alignment",   "higher",   95,      90,      "pct",
  "dedup_rate",           "Dedup rate",              "alignment",   "higher",   35,      25,      "pct",
  "total_counts",         "Total counts",            "libsize",     "higher",   5e6,     2.5e6,   "count",
  "top100_frac",          "Counts in top 100 genes", "complexity",  "lower",    40,      50,      "pct",
  "detected_genes",       "Detected genes",          "detection",   "higher",   22000,   19000,   "count",
  "genes_min10",          "Genes with >= 10 counts", "detection",   "higher",   13000,   10000,   "count",
  "protein_coding_pct",   "Protein-coding reads",    "biotype",     "higher",   85,      80,      "pct",
  "rRNA_pct",             "rRNA reads",              "biotype",     "lower",    5,       10,      "pct",
  "mito_pct",             "Mitochondrial reads",     "biotype",     "lower",    10,      12,      "pct",
  "within_condition_cor", "Within-condition r",      "clustering",  "higher",   0.95,    0.92,    "ratio",
  "pca_centroid_dist",    "PCA centroid distance",   "pca",         "lower",    5,       8,       "ratio",
  "median_abs_M",         "Replicate deviation",     "replicate",   "lower",    0.30,    0.45,    "ratio"
)

# Non-threshold QC knobs
QC_PARAMS <- list(
  filter_min_count = 10,   # filter_genes(): min counts ...
  filter_min_rep   = 3,    # ... in at least this many samples, before vst
  complexity_top_n = 100,  # genes summed for top100_frac
  pca_n_feats      = 500,  # top-variable genes for PCA, clustering, centroid distance
  pca_n_pcs        = 5,    # PCs in the scatter matrix
  detection_min    = 10    # count threshold for genes_min10
)

# ---- DESeqDataSet ----
# Counts rounded to match the limma fit. make_dds() is unused: its required
# ah_record argument downloads an EnsDb at runtime. plot_biotypes() needs only
# gene_name and gene_biotype in rowData.
cnt <- round(as.matrix(counts[, sample_map$sample, drop = FALSE]))

dds <- DESeq2::DESeqDataSetFromMatrix(
  cnt,
  colData = as.data.frame(sample_map),
  design  = ~1
)

bt_idx  <- match(rownames(dds), vendor_genes$gene_id)
biotype <- vendor_genes$gene_biotype[bt_idx]
symbol  <- vendor_genes$gene_name[bt_idx]

SummarizedExperiment::rowData(dds)$gene_name    <- symbol
SummarizedExperiment::rowData(dds)$gene_biotype <- biotype

vsd <- DESeq2::vst(
  RNAseqQC::filter_genes(dds,
                         min_count = QC_PARAMS$filter_min_count,
                         min_rep   = QC_PARAMS$filter_min_rep),
  blind = TRUE
)

# ---- Per-sample metrics ----
lib_size <- colSums(cnt)

# Percent of each library on the selected genes. NA, not 0, when the annotation
# carries none of the requested biotypes.
read_share <- function(keep) {
  keep <- !is.na(keep) & keep
  if (!any(keep)) return(rep(NA_real_, ncol(cnt)))
  colSums(cnt[keep, , drop = FALSE]) / lib_size * 100
}

qc_summary <- tibble::tibble(
  sample             = colnames(cnt),
  total_counts       = lib_size,
  detected_genes     = colSums(cnt > 0),
  genes_min10        = colSums(cnt >= QC_PARAMS$detection_min),
  top100_frac        = apply(cnt, 2, function(x) {
    sum(sort(x, decreasing = TRUE)[seq_len(QC_PARAMS$complexity_top_n)]) / sum(x) * 100
  }),
  rRNA_pct           = read_share(biotype %in% c("rRNA", "rRNA_pseudogene")),
  mito_pct           = read_share(grepl("^MT-", symbol)),
  protein_coding_pct = read_share(biotype == "protein_coding")
) |>
  dplyr::left_join(
    dplyr::select(sample_map, sample, line, stim, condition, experiment),
    by = "sample"
  )

# Vendor alignment stats, joined on plasmid code. Pre-normalisation, so a small
# library still shows.
if (file.exists(results_zip)) {
  vk <- vendor_key()

  map_stats <- utils::read.csv(
    unz(results_zip, "FCMSVB-mapping-stats-reads.csv"),
    check.names = FALSE
  ) |>
    dplyr::rename(plasmid = 1, uniquely_mapped = "Uniquely Mapped") |>
    dplyr::left_join(dplyr::select(sample_map, plasmid, sample), by = "plasmid") |>
    dplyr::select(sample, uniquely_mapped)

  qc_summary <- dplyr::left_join(qc_summary, map_stats, by = "sample")

  # Per-sample mapping TSVs, keyed by idx/plasmid as in 1_df_count.R. Two columns,
  # Metric and Value; each metric becomes a snake_case column.
  dup_stats <- do.call(rbind, lapply(seq_len(nrow(vk)), function(i) {
    m <- utils::read.delim(
      unz(results_zip,
          sprintf("FCMSVB_mapping-stats/FCMSVB_%d_%s.tsv", vk$idx[i], vk$plasmid[i])),
      check.names = FALSE, stringsAsFactors = FALSE
    )
    v <- setNames(as.list(m$Value),
                  tolower(gsub("[^A-Za-z0-9]+", "_", m$Metric)))

    tibble::as_tibble(c(list(plasmid = vk$plasmid[i]), v))
  })) |>
    # Aliases for the names QC_THRESHOLDS uses
    dplyr::rename(
      dedup_unique_reads = dedup_uniquely_mapped_reads,
      dedup_rate         = dedup_mapped_reads_rate
    ) |>
    dplyr::left_join(dplyr::select(sample_map, plasmid, sample), by = "plasmid") |>
    dplyr::select(-plasmid)

  qc_summary <- dplyr::left_join(qc_summary, dup_stats, by = "sample")
}

# ---- vst-derived metrics ----
# Correlation and replicate agreement measured within condition, against a
# sample's own replicates.
vmat <- SummarizedExperiment::assay(vsd)[, qc_summary$sample, drop = FALSE]
grp  <- as.character(qc_summary$condition)

cor_mat <- stats::cor(vmat, method = "pearson")
diag(cor_mat) <- NA

peers <- lapply(seq_along(grp), function(i) setdiff(which(grp == grp[i]), i))

qc_summary$within_condition_cor <- vapply(seq_along(peers), function(i) {
  if (length(peers[[i]]) == 0) NA_real_
  else stats::median(cor_mat[i, peers[[i]]], na.rm = TRUE)
}, numeric(1))

qc_summary$median_abs_M <- vapply(seq_along(peers), function(i) {
  p <- peers[[i]]
  if (length(p) == 0) return(NA_real_)
  ref <- if (length(p) == 1) vmat[, p] else matrixStats::rowMedians(vmat[, p, drop = FALSE])
  stats::median(abs(vmat[, i] - ref), na.rm = TRUE)
}, numeric(1))

# One PCA feeds both pca_centroid_dist and the PCA figure
top_var  <- order(matrixStats::rowVars(vmat), decreasing = TRUE)[seq_len(QC_PARAMS$pca_n_feats)]
pca      <- stats::prcomp(t(vmat[top_var, ]), scale. = FALSE)
pca_var  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc_x     <- pca$x[, 1:2]

centroid <- t(vapply(split(as.data.frame(pc_x), grp), colMeans, numeric(2)))[grp, , drop = FALSE]
qc_summary$pca_centroid_dist <- sqrt(rowSums((pc_x - centroid)^2))
qc_summary$PC1 <- pc_x[, 1]
qc_summary$PC2 <- pc_x[, 2]

# ---- Evaluate the checks ----
qc_checks <- dplyr::filter(QC_THRESHOLDS, metric %in% colnames(qc_summary))

qc_status <- do.call(rbind, lapply(seq_len(nrow(qc_checks)), function(i) {
  ck <- qc_checks[i, ]
  v  <- qc_summary[[ck$metric]]

  status <- if (ck$direction == "higher") {
    ifelse(v < ck$fail, "FAIL", ifelse(v < ck$warn, "WARN", "PASS"))
  } else {
    ifelse(v > ck$fail, "FAIL", ifelse(v > ck$warn, "WARN", "PASS"))
  }

  tibble::tibble(
    sample     = qc_summary$sample,
    experiment = qc_summary$experiment,
    metric     = ck$label,
    category   = ck$category,
    unit       = ck$unit,
    direction  = ck$direction,
    value      = v,
    warn_at    = ck$warn,
    fail_at    = ck$fail,
    status     = ifelse(is.na(v), "NA", status)
  )
})) |>
  dplyr::mutate(
    sample = factor(sample, levels = sample_map$sample),
    metric = factor(metric, levels = qc_checks$label),
    status = factor(status, levels = c("PASS", "WARN", "FAIL", "NA"))
  )

qc_fail_samples <- sort(unique(as.character(
  qc_status$sample[qc_status$status == "FAIL"]
)))

n_status <- table(qc_status$status)
message(sprintf("QC: %d pass, %d warn, %d fail over %d samples x %d checks",
                n_status[["PASS"]], n_status[["WARN"]], n_status[["FAIL"]],
                dplyr::n_distinct(qc_status$sample), nrow(qc_checks)))
if (length(qc_fail_samples) > 0) {
  message("FAIL on at least one check: ", paste(qc_fail_samples, collapse = ", "))
}

# ---- Plots ----
# bg = "white" throughout; the ggsave default is transparent.
fig_qc <- ensure_dir(dir_fig, "qc")

png_out <- function(p, name, width = 9, height = 6) {
  ggsave(file.path(fig_qc, name), p, width = width, height = height,
         dpi = 300, bg = "white")
}

# Sample-labelled x axis, shared by the bar and PCA figures
x_samples <- theme(
  axis.text.x        = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                    colour = "black", size = 6),
  panel.grid.major.x = element_blank()
)

# Per-sample bars, one facet per check. Threshold lines drawn only where they fall
# inside the observed range; an out-of-range cutoff flattens the panel.
guides_df <- qc_status |>
  dplyr::group_by(metric) |>
  dplyr::summarise(lo = min(value, na.rm = TRUE), hi = max(value, na.rm = TRUE),
                   warn_at = warn_at[1], fail_at = fail_at[1], .groups = "drop") |>
  tidyr::pivot_longer(c(warn_at, fail_at), names_to = "level", values_to = "y") |>
  dplyr::filter(y >= lo, y <= hi)

qc_bars <- ggplot(qc_status, aes(x = sample, y = value, fill = status)) +
  geom_col() +
  geom_hline(data = guides_df, aes(yintercept = y, linetype = level),
             colour = "grey30", linewidth = 0.4) +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = qc_status_cols, drop = FALSE) +
  scale_linetype_manual(values = c(warn_at = "dashed", fail_at = "dotted"),
                        labels = c(warn_at = "warn", fail_at = "fail")) +
  scale_y_continuous(labels = scales::label_number(scale_cut = scales::cut_short_scale())) +
  labs(x = NULL, y = NULL, fill = NULL, linetype = NULL) +
  theme_bw() +
  x_samples +
  theme(legend.position = "top")

png_out(qc_bars, "qc_metrics_by_sample.png", 15, 18)

# PCA from the same prcomp as pca_centroid_dist. Only flagged samples are labelled.
pca_df <- qc_summary |>
  dplyr::mutate(flagged = sample %in% qc_status$sample[qc_status$status != "PASS"])

qc_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = line, shape = stim)) +
  geom_point(size = 3, colour = "grey20", alpha = 0.9) +
  ggrepel::geom_text_repel(
    data = dplyr::filter(pca_df, flagged),
    aes(x = PC1, y = PC2, label = sample), size = 2.6, colour = "grey15",
    min.segment.length = 0, max.overlaps = Inf, inherit.aes = FALSE,
    box.padding = 0.4
  ) +
  scale_fill_manual(values = line_cols) +
  scale_shape_manual(values = stim_shapes) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"),
       fill = NULL, shape = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  theme_bw()

png_out(qc_pca, "qc_pca.png", 8, 6)

# RNAseqQC distribution views, run-level and unlabelled
png_out(RNAseqQC::plot_total_counts(dds),                    "qc_total_counts.png")
png_out(RNAseqQC::plot_library_complexity(dds, show_progress = FALSE),
        "qc_library_complexity.png")
png_out(RNAseqQC::plot_gene_detection(dds),                  "qc_gene_detection.png")
png_out(RNAseqQC::plot_biotypes(dds),                        "qc_biotypes.png")
png_out(RNAseqQC::mean_sd_plot(vsd),                         "qc_mean_sd.png", 7, 5)
png_out(RNAseqQC::plot_pca_scatters(vsd, n_PCs = QC_PARAMS$pca_n_pcs,
                                    n_feats = QC_PARAMS$pca_n_feats,
                                    color_by = "line", shape_by = "stim"),
        "qc_pca_scatters.png", 11, 10)

# ComplexHeatmap, not ggplot: needs its own device
png(file.path(fig_qc, "qc_sample_clustering.png"),
    width = 10, height = 9, units = "in", res = 300, bg = "white")
ComplexHeatmap::draw(
  RNAseqQC::plot_sample_clustering(vsd,
                                   n_feats   = QC_PARAMS$pca_n_feats,
                                   anno_vars = c("line", "stim"),
                                   distance  = "euclidean")
)
dev.off()

# Status grid: one tile per sample x check, value in the tile, bold when flagged.
# Metrics reversed to match QC_THRESHOLDS order top to bottom.
qc_tile <- ggplot(qc_status, aes(x = sample, y = metric, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.7) +
  geom_text(
    aes(label = fmt_qc(value, unit, compact = TRUE),
        fontface = ifelse(status %in% c("WARN", "FAIL"), "bold", "plain")),
    size = 2.3, colour = "grey15"
  ) +
  facet_grid(. ~ experiment, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = qc_status_cols, drop = FALSE) +
  scale_y_discrete(limits = rev(levels(qc_status$metric))) +
  labs(x = NULL, y = NULL, fill = NULL) +
  theme_bw() +
  theme(
    axis.text.x     = element_text(angle = 90, hjust = 1, vjust = 0.5, colour = "black"),
    axis.text.y     = element_text(colour = "black"),
    panel.grid      = element_blank(),
    legend.position = "top"
  )

ggsave(file.path(fig_qc, "qc_status_grid.png"),
       qc_tile, width = 14, height = 1 + 0.42 * nrow(qc_checks), dpi = 300, bg = "white")

# ---- Report ----
write_qc_report(
  path              = file.path(ensure_dir(dir_res), "qc_report.md"),
  qc_summary        = qc_summary,
  qc_status         = qc_status,
  thresholds        = QC_THRESHOLDS,
  params            = QC_PARAMS,
  prov              = capture_provenance(),
  counts_checkpoint = basename(tail(list.files(dir_rds,
                                               pattern = "^annotated_counts_v\\d+\\.rds$",
                                               full.names = TRUE), 1))
)
