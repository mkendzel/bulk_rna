# Cell-composition / maturity figures. Scores marker modules per sample and
# tests C9 against the pooled controls. Writes results/composition_module_stats.tsv
# and the module score checkpoints 9_composition_adjusted.R reads.
# Run from the repo root, after 3_DE_limma.R.
#
# z is taken within experiment, so spinal and cort scores are not on a shared
# scale.

# ---- Library ----
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(patchwork)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Data load ----
tt          <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("tt_", e, "_annotated"), dir = dir_rds)), EXPERIMENTS)
logcpm      <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("logcpm_", e), dir = dir_rds)), EXPERIMENTS)
sample_meta <- setNames(lapply(EXPERIMENTS, function(e) load_checkpoint(paste0("sample_meta_", e), dir = dir_rds)), EXPERIMENTS)

# =============================================================================
# Marker modules
# =============================================================================
# Four modules spanning the differentiation axis of a spinal/cortical organoid.
# Non-overlapping by construction - the stopifnot below enforces it. Progenitor
# is radial glia (VIM, FABP7, SLC1A3, HES1/5), not generic stem markers.
MATURITY_MODULES <- list(
  Proliferation = c("MKI67", "TOP2A", "CCNB1", "PCNA", "CDK1", "BIRC5", "CENPF", "AURKB"),
  Progenitor    = c("PAX6", "SOX2", "NES", "VIM", "HES1", "HES5", "FABP7", "SLC1A3"),
  Neuron        = c("MAP2", "RBFOX3", "SYP", "STMN2", "DCX", "TUBB3", "SNAP25", "NEFL", "ELAVL3", "SYT1"),
  Astrocyte     = c("AQP4", "GFAP", "S100B", "ALDH1L1", "SLC1A2", "AGT")
)

MATURE_MODULES   <- c("Neuron", "Astrocyte")
IMMATURE_MODULES <- c("Proliferation", "Progenitor")

stopifnot(!any(duplicated(unlist(MATURITY_MODULES))))

CASE_LINE  <- "C9"
CTRL_LINES <- c("Ctrl3", "Ctrl2")

dir_comp <- file.path("composition")

# =============================================================================
# Stats
# =============================================================================
# C9 minus pooled-control mean per module, signed like the genotype contrasts in
# 0_config.R (positive = higher in C9). p is a nominal, unadjusted t-test over
# organoid batches; safe_t_p() is in R/plot_guards.R.
case_vs_control <- function(d, scope) {
  d |>
    dplyr::mutate(grp = dplyr::case_when(
      line == CASE_LINE   ~ "case",
      line %in% CTRL_LINES ~ "control",
      TRUE                ~ NA_character_
    )) |>
    dplyr::filter(!is.na(grp)) |>
    dplyr::group_by(module) |>
    dplyr::summarise(
      scope    = scope,
      n_case   = sum(grp == "case"),
      n_ctrl   = sum(grp == "control"),
      delta    = mean(score[grp == "case"], na.rm = TRUE) -
                 mean(score[grp == "control"], na.rm = TRUE),
      p_nominal = safe_t_p(score[grp == "case"], score[grp == "control"]),
      .groups  = "drop"
    ) |>
    dplyr::relocate(scope, .after = module)
}

# =============================================================================
# Figures
# =============================================================================
# Shared theme for the composition panels: line on fill, stimulus on shape
comp_theme <- function() {
  theme_bw() +
    theme(
      panel.grid.major.x = element_blank(),
      strip.background   = element_rect(fill = "grey92", colour = NA),
      strip.text         = element_text(face = "bold", size = 9),
      legend.position    = "bottom",
      legend.box         = "horizontal"
    )
}

# Panel A: one facet per module, boxplot per line with every sample overlaid
module_panel <- function(scores, n_lab) {

  if (nrow(scores) == 0L) return(empty_plot("Module scores", note = "No markers detected"))

  ggplot(scores, aes(line, score)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_boxplot(aes(fill = line), width = 0.6, outlier.shape = NA,
                 alpha = 0.55, linewidth = 0.3) +
    geom_point(aes(shape = stim), position = position_jitter(width = 0.14, height = 0),
               size = 1.9, fill = "white", stroke = 0.4) +
    facet_wrap(~ module, nrow = 1, labeller = labeller(module = n_lab)) +
    scale_fill_manual(values = line_cols, guide = "none") +
    scale_shape_manual(values = stim_shapes, name = NULL) +
    labs(x = NULL, y = "Module score (mean z)",
         title = "A. Marker module scores, one point per sample") +
    comp_theme()
}

# Panel B: the four modules collapsed to one maturity index, samples sorted
maturity_panel <- function(wide) {

  if (nrow(wide) == 0L) return(empty_plot("Maturity index", note = "No markers detected"))

  d <- dplyr::mutate(wide, sample = factor(sample, levels = sample[order(maturity)]))

  ggplot(d, aes(sample, maturity)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_segment(aes(xend = sample, y = 0, yend = maturity, colour = line),
                 linewidth = 0.5) +
    geom_point(aes(fill = line, shape = stim), size = 2.6, stroke = 0.4) +
    scale_fill_manual(values = line_cols, name = NULL) +
    scale_colour_manual(values = line_cols, guide = "none") +
    scale_shape_manual(values = stim_shapes, name = NULL) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    labs(x = NULL,
         y = "Maturity index\n(neuron + astro) - (prolif + progenitor)",
         title = "B. One maturity axis per sample, sorted") +
    comp_theme() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6))
}

# Figure 2: per-gene z-scores behind the module scores, one tile per gene x sample
marker_heatmap <- function(z_long, title, subtitle) {

  if (nrow(z_long) == 0L) return(empty_plot(title, subtitle, "No markers detected"))

  ord <- z_long |>
    dplyr::distinct(sample, line, stim) |>
    dplyr::arrange(line, stim, sample) |>
    dplyr::pull(sample)

  # Genes sorted within module by C9-minus-control difference
  gene_ord <- z_long |>
    dplyr::group_by(module, SYMBOL) |>
    dplyr::summarise(
      d = mean(z[line == CASE_LINE], na.rm = TRUE) -
          mean(z[line %in% CTRL_LINES], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(module, d) |>
    dplyr::pull(SYMBOL)

  d <- z_long |>
    dplyr::mutate(
      sample = factor(sample, levels = ord),
      SYMBOL = factor(SYMBOL, levels = gene_ord),
      z      = pmax(pmin(z, 2), -2)   # clip so one outlier does not eat the ramp
    )

  ggplot(d, aes(sample, SYMBOL, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    facet_grid(module ~ line, scales = "free", space = "free") +
    scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                         midpoint = 0, limits = c(-2, 2),
                         name = "z (clipped)") +
    labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    comp_theme() +
    theme(
      axis.text.x     = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 6),
      axis.text.y     = element_text(size = 7),
      panel.grid      = element_blank(),
      legend.position = "right"
    )
}

# =============================================================================
# Build
# =============================================================================
stats_all <- list()

for (e in EXPERIMENTS) {

  lc  <- logcpm[[e]]
  sm  <- sample_meta[[e]]
  map <- symbol_map(tt[[e]])

  z_long <- module_expr_z(lc, sm, map, MATURITY_MODULES)

  lost <- setdiff(unlist(MATURITY_MODULES), unique(z_long$SYMBOL))
  if (length(lost) > 0) {
    message(e, " / maturity markers not detected: ", paste(lost, collapse = ", "))
  }

  z_long <- dplyr::mutate(z_long, line = droplevels(line), stim = droplevels(stim))

  scores <- module_score(z_long)
  wide   <- module_score_wide(scores, MATURE_MODULES, IMMATURE_MODULES)

  save_checkpoint(scores, paste0("module_scores_", e), dir = dir_rds)
  save_checkpoint(wide,   paste0("module_scores_wide_", e), dir = dir_rds)

  # Facet strips carry the detected / total marker count
  n_by_mod <- scores |>
    dplyr::distinct(module, n_genes) |>
    dplyr::mutate(want = lengths(MATURITY_MODULES)[as.character(module)])
  n_lab <- setNames(
    sprintf("%s\n%d/%d genes", n_by_mod$module, n_by_mod$n_genes, n_by_mod$want),
    as.character(n_by_mod$module)
  )

  has_case <- CASE_LINE %in% levels(scores$line)

  if (has_case) {
    st <- dplyr::bind_rows(
      case_vs_control(scores, "all_samples"),
      case_vs_control(dplyr::filter(scores, stim == "Mock"), "mock_only")
    ) |>
      dplyr::mutate(experiment = e, .before = 1)
    stats_all[[e]] <- st
    print(st)

    hit <- st |> dplyr::filter(scope == "all_samples")
    sub_a <- paste0(
      "C9 - control (mean z): ",
      paste(sprintf("%s %+.2f", hit$module, hit$delta), collapse = "; ")
    )
  } else {
    sub_a <- "Single line in this experiment - no genotype comparison"
  }

  cap <- paste(
    "Scores are the mean z-score of each module's detected markers, z taken across",
    "all samples of this experiment. A module score tracks the fraction of a cell",
    "type but is not a cell count, and replicates are organoid batches from one",
    "line per genotype - so this measures the composition difference, it does not",
    "attribute it to C9orf72.",
    sep = "\n"
  )

  fig <- (module_panel(scores, n_lab) / maturity_panel(wide)) +
    plot_layout(heights = c(1, 1.15)) +
    plot_annotation(
      title    = paste0("Organoid composition and maturity - ", e),
      subtitle = sub_a,
      caption  = cap,
      theme    = theme(
        plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey30", size = 9),
        plot.caption  = element_text(colour = "grey40", size = 7, hjust = 0)
      )
    )

  save_fig(fig, paste0("composition_modules_", e),
           width = 11, height = 9, subdir = dir_comp)

  save_fig(
    marker_heatmap(z_long,
                   paste0("Maturity markers by sample - ", e),
                   "z within experiment, clipped at +/-2; genes sorted within module by C9 - control"),
    paste0("composition_markers_", e),
    width = max(8, 0.42 * dplyr::n_distinct(z_long$sample) + 5),
    height = 9, subdir = dir_comp)
}

if (length(stats_all) > 0) {
  out <- file.path(ensure_dir(dir_res), "composition_module_stats.tsv")
  utils::write.table(dplyr::bind_rows(stats_all), out,
                     sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", out)
}

message("8_composition.R complete")
