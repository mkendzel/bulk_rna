# Shared config for WNV_ALS_R01_2026. Sourced by scripts 1-7.
# Run from the repo root.

library(dplyr)
library(tibble)
library(tidyr)

# ---- Paths ----
PROJ    <- "analysis/WNV_ALS_R01_2026"
dir_raw <- file.path(PROJ, "data", "raw")
dir_rds <- file.path(PROJ, "data", "r_objects")
dir_fig <- file.path(PROJ, "figures")
dir_res <- file.path(PROJ, "results")

expr_path   <- file.path(dir_raw, "FCMSVB-expression-matrix.tsv")
results_zip <- file.path(dir_raw, "FCMSVB_results.zip")

# Create a directory and return its path
ensure_dir <- function(...) {
  p <- file.path(...)
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

# ---- Significance cutoffs ----
P_CUT   <- 0.05
LFC_CUT <- 1.5

# ---- Sample map ----
# Plasmid code -> real sample name, from notes/SM072026_bulk RNAseq.xlsx.
# Code grammar: <LINE>-<STIM><rep>_<globalIndex>
#   line  C9 | CT3 -> Ctrl3 | CT2 -> Ctrl2 | CO -> Cort
#   stim  M -> Mock | b -> IFNb | g -> IFNg | TX -> WNV
sample_map <- tibble::tribble(
  ~plasmid,     ~sample,        ~line,   ~stim,  ~rep, ~experiment,
  # C9_1 ALS spinal organoids
  "C9-M1_1",    "C9_Mock_1",    "C9",    "Mock",   1L, "spinal",
  "C9-M2_2",    "C9_Mock_2",    "C9",    "Mock",   2L, "spinal",
  "C9-M3_3",    "C9_Mock_3",    "C9",    "Mock",   3L, "spinal",
  "C9-b1_4",    "C9_IFNb_1",    "C9",    "IFNb",   1L, "spinal",
  "C9-b2_5",    "C9_IFNb_2",    "C9",    "IFNb",   2L, "spinal",
  "C9-b3_6",    "C9_IFNb_3",    "C9",    "IFNb",   3L, "spinal",
  "C9-g1_7",    "C9_IFNg_1",    "C9",    "IFNg",   1L, "spinal",
  "C9-g2_8",    "C9_IFNg_2",    "C9",    "IFNg",   2L, "spinal",
  "C9-g3_9",    "C9_IFNg_3",    "C9",    "IFNg",   3L, "spinal",
  "C9-TX1_10",  "C9_WNV_1",     "C9",    "WNV",    1L, "spinal",
  "C9-TX2_11",  "C9_WNV_2",     "C9",    "WNV",    2L, "spinal",
  "C9-TX3_12",  "C9_WNV_3",     "C9",    "WNV",    3L, "spinal",
  # Ctrl_3 healthy control. No Mock 3.
  "CT3-M1_13",  "Ctrl3_Mock_1", "Ctrl3", "Mock",   1L, "spinal",
  "CT3-M2_14",  "Ctrl3_Mock_2", "Ctrl3", "Mock",   2L, "spinal",
  "CT3-b1_15",  "Ctrl3_IFNb_1", "Ctrl3", "IFNb",   1L, "spinal",
  "CT3-b2_16",  "Ctrl3_IFNb_2", "Ctrl3", "IFNb",   2L, "spinal",
  "CT3-b3_17",  "Ctrl3_IFNb_3", "Ctrl3", "IFNb",   3L, "spinal",
  "CT3-g1_18",  "Ctrl3_IFNg_1", "Ctrl3", "IFNg",   1L, "spinal",
  "CT3-g2_19",  "Ctrl3_IFNg_2", "Ctrl3", "IFNg",   2L, "spinal",
  "CT3-g3_20",  "Ctrl3_IFNg_3", "Ctrl3", "IFNg",   3L, "spinal",
  "CT3-TX1_21", "Ctrl3_WNV_1",  "Ctrl3", "WNV",    1L, "spinal",
  "CT3-TX2_22", "Ctrl3_WNV_2",  "Ctrl3", "WNV",    2L, "spinal",
  "CT3-TX3_23", "Ctrl3_WNV_3",  "Ctrl3", "WNV",    3L, "spinal",
  # Ctrl_2 healthy control. Single IFNb and IFNg replicate.
  "CT2-M1_24",  "Ctrl2_Mock_1", "Ctrl2", "Mock",   1L, "spinal",
  "CT2-M2_25",  "Ctrl2_Mock_2", "Ctrl2", "Mock",   2L, "spinal",
  "CT2-M3_26",  "Ctrl2_Mock_3", "Ctrl2", "Mock",   3L, "spinal",
  "CT2-b1_27",  "Ctrl2_IFNb_1", "Ctrl2", "IFNb",   1L, "spinal",
  "CT2-g1_28",  "Ctrl2_IFNg_1", "Ctrl2", "IFNg",   1L, "spinal",
  "CT2-TX1_29", "Ctrl2_WNV_1",  "Ctrl2", "WNV",    1L, "spinal",
  "CT2-TX2_30", "Ctrl2_WNV_2",  "Ctrl2", "WNV",    2L, "spinal",
  "CT2-TX3_31", "Ctrl2_WNV_3",  "Ctrl2", "WNV",    3L, "spinal",
  # Cortical organoids, uninfected. Separate experiment, no WNV arm.
  "CO-M1_32",   "Cort_Mock_1",  "Cort",  "Mock",   1L, "cort",
  "CO-M2_33",   "Cort_Mock_2",  "Cort",  "Mock",   2L, "cort",
  "CO-M3_34",   "Cort_Mock_3",  "Cort",  "Mock",   3L, "cort",
  "CO-b1_35",   "Cort_IFNb_1",  "Cort",  "IFNb",   1L, "cort",
  "CO-b2_36",   "Cort_IFNb_2",  "Cort",  "IFNb",   2L, "cort",
  "CO-b3_37",   "Cort_IFNb_3",  "Cort",  "IFNb",   3L, "cort",
  "CO-g1_38",   "Cort_IFNg_1",  "Cort",  "IFNg",   1L, "cort",
  "CO-g2_39",   "Cort_IFNg_2",  "Cort",  "IFNg",   2L, "cort",
  "CO-g3_40",   "Cort_IFNg_3",  "Cort",  "IFNg",   3L, "cort"
) |>
  dplyr::mutate(
    condition = paste(line, stim, sep = "_"),
    idx       = as.integer(sub("^.*_", "", plasmid))  # fallback key for column matching
  )

stopifnot(
  nrow(sample_map) == 40L,
  !any(duplicated(sample_map$plasmid)),
  !any(duplicated(sample_map$sample)),
  identical(sort(sample_map$idx), 1:40)
)

# Vendor column index -> plasmid code, parsed from the per-sample filenames in the
# results zip (FCMSVB_mapping-stats/FCMSVB_<idx>_<plasmid>.tsv).
# Returns NULL when the zip is absent.
vendor_key <- function(zip = results_zip) {
  if (!file.exists(zip)) return(NULL)
  nm <- utils::unzip(zip, list = TRUE)$Name
  m  <- regmatches(nm, regexec("mapping-stats/FCMSVB_(\\d+)_(.+)\\.tsv$", nm))
  m  <- m[lengths(m) == 3L]
  if (length(m) == 0L) return(NULL)
  tibble::tibble(
    idx     = as.integer(vapply(m, `[`, character(1), 2L)),
    plasmid = vapply(m, `[`, character(1), 3L)
  ) |>
    dplyr::arrange(idx)
}

# ---- Factor levels ----
LINE_LEVELS <- c("C9", "Ctrl3", "Ctrl2", "Cort")
STIM_LEVELS <- c("Mock", "IFNb", "IFNg", "WNV")
EXPERIMENTS <- c("spinal", "cort")

sample_map$line <- factor(sample_map$line, levels = LINE_LEVELS)
sample_map$stim <- factor(sample_map$stim, levels = STIM_LEVELS)

# Ordered line then stim; fixes design matrix column order
cond_levels <- sample_map |>
  dplyr::distinct(line, stim, condition) |>
  dplyr::arrange(line, stim) |>
  dplyr::pull(condition)

sample_map$condition <- factor(sample_map$condition, levels = cond_levels)

# Replicates per group. Ctrl2_IFNb and Ctrl2_IFNg are n = 1.
cond_n <- sample_map |>
  dplyr::count(condition, name = "n") |>
  tibble::deframe()

# ---- Contrast registry ----
# One row per contrast. `expr` goes to limma::makeContrasts; scripts 3-7 join on
# `name` for grouping metadata and plot labels.

# Each stimulus vs its own line's Mock
.within <- dplyr::bind_rows(
  tidyr::expand_grid(
    line = c("C9", "Ctrl3", "Ctrl2"),
    stim = c("IFNb", "IFNg", "WNV")
  ) |>
    dplyr::mutate(experiment = "spinal"),
  tidyr::expand_grid(
    line = "Cort",
    stim = c("IFNb", "IFNg")
  ) |>
    dplyr::mutate(experiment = "cort")
) |>
  dplyr::transmute(
    name       = sprintf("%s_%s_vs_Mock", line, stim),
    type       = "within_line",
    experiment = experiment,
    line       = line,
    ref_line   = NA_character_,
    stim       = stim,
    expr       = sprintf("%s_%s - %s_Mock", line, stim, line),
    label      = sprintf("%s: %s vs Mock", line, stim)
  )

# ALS line vs each healthy control line at matched stimulus
.genotype <- tidyr::expand_grid(
  ref_line = c("Ctrl3", "Ctrl2"),
  stim     = STIM_LEVELS
) |>
  dplyr::transmute(
    name       = sprintf("C9_vs_%s_%s", ref_line, stim),
    type       = "genotype",
    experiment = "spinal",
    line       = "C9",
    ref_line   = ref_line,
    stim       = stim,
    expr       = sprintf("C9_%s - %s_%s", stim, ref_line, stim),
    label      = sprintf("%s: C9 vs %s", stim, ref_line)
  )

# Difference of differences: ALS response to stimulus vs control response
.interaction <- tidyr::expand_grid(
  ref_line = c("Ctrl3", "Ctrl2"),
  stim     = c("IFNb", "IFNg", "WNV")
) |>
  dplyr::transmute(
    name       = sprintf("%s_response_C9_vs_%s", stim, ref_line),
    type       = "interaction",
    experiment = "spinal",
    line       = "C9",
    ref_line   = ref_line,
    stim       = stim,
    expr       = sprintf("(C9_%s - C9_Mock) - (%s_%s - %s_Mock)",
                         stim, ref_line, stim, ref_line),
    label      = sprintf("%s response: C9 vs %s", stim, ref_line)
  )

contrast_registry <- dplyr::bind_rows(.within, .genotype, .interaction)

# Smallest group size in each contrast; used to flag unreplicated results
.groups_in_expr <- function(e) {
  toks <- unlist(strsplit(e, "[^A-Za-z0-9_]+"))
  intersect(toks, cond_levels)
}

contrast_registry$min_n <- vapply(
  contrast_registry$expr,
  function(e) min(cond_n[.groups_in_expr(e)]),
  integer(1)
)

contrast_registry <- contrast_registry |>
  dplyr::mutate(
    type = factor(type, levels = c("within_line", "genotype", "interaction")),
    line = factor(line, levels = LINE_LEVELS),
    stim = factor(stim, levels = STIM_LEVELS)
  ) |>
  dplyr::arrange(experiment, type, line, ref_line, stim)

rm(.within, .genotype, .interaction)

stopifnot(
  nrow(contrast_registry) == 25L,
  !any(duplicated(contrast_registry$name)),
  sum(contrast_registry$experiment == "spinal") == 23L,
  sum(contrast_registry$experiment == "cort") == 2L
)

# Contrast names for one experiment, in registry order
contrasts_for <- function(experiment) {
  contrast_registry$name[contrast_registry$experiment == experiment]
}

# Contrast labels for one experiment, in registry order. Shared x-axis level set
# for the figures in scripts 5-7. The `registry` default resolves lazily against
# globalenv(), so a script that overwrites contrast_registry with the checkpoint
# gets the min_n 3_DE_limma.R recomputed after drop_samples.
contrast_labels <- function(experiment, registry = contrast_registry) {
  exp_name <- experiment
  registry |>
    dplyr::filter(experiment == exp_name) |>
    dplyr::arrange(type, line, ref_line, stim) |>
    dplyr::pull(label) |>
    unique()
}

# Registry row for a contrast name; an all-NA row for an unknown name.
reg_of <- function(cn, registry = contrast_registry) {
  registry[match(cn, registry$name), ]
}

# ---- Palettes ----
line_cols <- c(C9 = "#D62728", Ctrl3 = "#1F77B4", Ctrl2 = "#17BECF", Cort = "#7F7F7F")

# Tint a colour toward white
make_shades <- function(hex, n = 4) {
  grDevices::colorRampPalette(c("#FFFFFF", hex))(n + 1)[-1]
}

# One colour per condition: line hue, stimulus intensity
cond_cols <- unlist(lapply(LINE_LEVELS, function(l) {
  stims <- levels(droplevels(sample_map$stim[sample_map$line == l]))
  setNames(make_shades(line_cols[[l]], length(stims)),
           paste(l, stims, sep = "_"))
}))
cond_cols <- cond_cols[cond_levels]

stim_shapes <- c(Mock = 21, IFNb = 22, IFNg = 23, WNV = 24)

type_cols <- c(within_line = "#4C72B0", genotype = "#DD8452", interaction = "#55A868")

# QC status tiles in script 2
qc_status_cols <- c(PASS = "#EFF3EF", WARN = "#FBE0A6", FAIL = "#F2B4AC", "NA" = "#E6E6E6")
