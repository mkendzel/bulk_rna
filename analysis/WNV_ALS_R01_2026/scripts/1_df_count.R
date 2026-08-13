# Import the Plasmidsaurus matrix, rename columns to real sample names, checkpoint.
# Nothing is filtered or dropped.
#
# Run once per delivery. Scripts 2 and 3 load the checkpoints instead of re-importing.
# Run from the repo root.

# ---- Libraries ----
library(dplyr)
library(tibble)

# ---- Load helper functions + project config ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ----- Import data ----
# Plasmidsaurus matrix: gene_id / gene_name / gene_biotype, then a <sample>_cpm and
# <sample>_count column per sample. Counts are fractional; rounded downstream.
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
# Vendor ships columns as FCMSVB_<idx>. Match on the vendor idx -> code map first,
# then the exact plasmid code, then the trailing index. The vendor map is checked
# against sample_map first, so a renumbered run errors instead of permuting labels.
vk <- vendor_key()

rename_by_vendor <- character(0)
if (!is.null(vk)) {
  chk <- dplyr::inner_join(vk, sample_map[, c("idx", "plasmid", "sample")],
                           by = "idx", suffix = c("_vendor", "_map"))
  stopifnot(
    nrow(chk) == nrow(sample_map),
    identical(chk$plasmid_vendor, chk$plasmid_map)
  )
  rename_by_vendor <- setNames(chk$sample, paste0("FCMSVB_", chk$idx))
}

rename_by_plasmid <- setNames(sample_map$sample, sample_map$plasmid)
rename_by_idx     <- setNames(sample_map$sample, as.character(sample_map$idx))

new_names <- vapply(colnames(counts), function(cn) {
  if (cn %in% names(rename_by_vendor))  return(rename_by_vendor[[cn]])
  if (cn %in% names(rename_by_plasmid)) return(rename_by_plasmid[[cn]])
  idx <- sub("^.*_", "", cn)
  if (grepl("^[0-9]+$", idx) && idx %in% names(rename_by_idx)) return(rename_by_idx[[idx]])
  NA_character_
}, character(1), USE.NAMES = FALSE)

# List failures before stopping
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

# ---- Checkpoints ----
save_checkpoint(counts, "annotated_counts", dir = dir_rds)
save_checkpoint(vendor_genes, "vendor_genes", dir = dir_rds)
