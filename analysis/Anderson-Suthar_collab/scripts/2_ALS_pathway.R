library(limma)
library(msigdbr)
library(dplyr)

invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# ── Load checkpoints ──
res_d30       <- load_checkpoint("res_d30_ALS_vs_Ctrl",     dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d50       <- load_checkpoint("res_d50_ALS_vs_Ctrl",     dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d75       <- load_checkpoint("res_d75_ALS_vs_Ctrl",     dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d120      <- load_checkpoint("res_d120_ALS_vs_Ctrl",    dir = "projects/Anderson-Suthar_collab/data/r_objects")
v             <- load_checkpoint("voom_hSpS_Ctrx_batch",    dir = "projects/Anderson-Suthar_collab/data/r_objects")
meta_filtered <- load_checkpoint("meta_filtered_hSpS_Ctrx", dir = "projects/Anderson-Suthar_collab/data/r_objects")

# ── WP2447 (ALS) gene set ──
wp2447_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") |>
  dplyr::filter(gs_exact_source == "WP2447")
als_syms <- unique(wp2447_df$gene_symbol)

# ── camera: correlation-aware competitive test per timepoint ──
design <- model.matrix(~ Sex + Batch + Stage * Line, data = meta_filtered)
colnames(design) <- make.names(colnames(design))

con <- makeContrasts(
  d30_ALS  = LineALS,
  d50_ALS  = LineALS + Staged50.LineALS,
  d75_ALS  = LineALS + Staged75.LineALS,
  d120_ALS = LineALS + Staged120.LineALS,
  levels   = design
)

als_idx <- which(rownames(v) %in% als_syms)

camera_res <- do.call(rbind, lapply(colnames(con), function(coef) {
  cam <- camera(v, index = als_idx, design = design, contrast = con[, coef],
                inter.gene.cor = NA)
  data.frame(timepoint = sub("_ALS$", "", coef), cam, row.names = NULL)
}))

# ── Magnitude and directionality summary from topTable results ──
res_list <- list(d30 = res_d30, d50 = res_d50, d75 = res_d75, d120 = res_d120)

summarise_set <- function(res, label) {
  idx      <- which(rownames(res) %in% als_syms)
  t_vals   <- res[["t"]][idx]
  lfc_vals <- res[["logFC"]][idx]

  data.frame(
    timepoint    = label,
    rms_t        = sqrt(mean(t_vals^2)),
    mean_abs_lfc = mean(abs(lfc_vals)),
    mean_t_net   = mean(t_vals),
    mean_lfc_net = mean(lfc_vals),
    prop_up      = mean(lfc_vals > 0),
    prop_down    = mean(lfc_vals < 0)
  )
}

set_summary <- do.call(rbind, Map(summarise_set, res_list, names(res_list)))

# ── Combine camera p-values with magnitude summary ──
results <- merge(camera_res, set_summary, by = "timepoint")
print(results)


# ---- Fry -----
# ── WP2447 (ALS) mixed (non-directional) set test ──
# uses objects from the script above: v, design, con, als_syms, res_d30/_d50/_d75/_d120

# config: set both only if voom_hSpS_Ctrx_batch was fit with duplicateCorrelation
dup_block <- NULL   # repeated-measures block variable; NULL = no correction
dup_cor   <- NULL   # consensus correlation from duplicateCorrelation()

als_idx  <- which(rownames(v) %in% als_syms)
res_list <- list(d30 = res_d30, d50 = res_d50, d75 = res_d75, d120 = res_d120)

# fry: self-contained rotation test; PValue.Mixed is the direction-agnostic result
fry_mixed <- function(coef) {
  args <- list(y = v, index = als_idx, design = design, contrast = con[, coef])
  if (!is.null(dup_block)) {
    args$block       <- dup_block
    args$correlation <- dup_cor
  }
  fr <- do.call(fry, args)
  data.frame(
    timepoint     = sub("_ALS$", "", coef),
    NGenes        = fr$NGenes,
    p_mixed       = fr$PValue.Mixed,
    p_directional = fr$PValue,
    direction     = fr$Direction
  )
}

fry_res <- do.call(rbind, lapply(colnames(con), fry_mixed))

# magnitude and split from the per-gene topTables (sign-agnostic, no cancellation)
summarise_set <- function(res, label) {
  idx      <- which(rownames(res) %in% als_syms)
  t_vals   <- res[["t"]][idx]
  lfc_vals <- res[["logFC"]][idx]
  data.frame(
    timepoint    = label,
    rms_t        = sqrt(mean(t_vals^2)),
    mean_abs_lfc = mean(abs(lfc_vals)),
    prop_up      = mean(lfc_vals > 0),
    prop_down    = mean(lfc_vals < 0)
  )
}

set_summary <- do.call(rbind, Map(summarise_set, res_list, names(res_list)))

results <- merge(fry_res, set_summary, by = "timepoint")
print(results)
