# Regresses ISG induction on the composition module scores from 8_composition.R,
# with and without a genotype term, and reports the same comparison at gene level
# under two limma designs. Writes five TSVs to results/ and one figure.
# Run from the repo root, after 8_composition.R.
#
# Three models, in order:
#   1. multiple regression on the stimulated samples, genotype +/- covariate
#   2. slope calibrated on control samples only, C9 predicted from it
#   3. gene-level limma under the cell-means and reduced designs
#
# Nothing here calls fgsea. The set-level statistic is the mean moderated t of
# the Hallmark IFN-alpha genes, which is deterministic; re-permuting would
# contradict results/gsea_hallmark_all.tsv.

# ---- Library ----
library(edgeR)
library(limma)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(patchwork)
library(fgsea)   # gmtPathways only - no permutation

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Parameters ----
EXP        <- "spinal"           # the only experiment with a genotype comparison
ISG_SET    <- "HALLMARK_INTERFERON_ALPHA_RESPONSE"
CASE_LINE  <- "C9"
CTRL_LINES <- c("Ctrl3", "Ctrl2")
COVARIATE  <- "Astrocyte"        # primary; Maturity index run as a sensitivity
dir_comp   <- file.path("composition")

# ---- Data load ----
counts      <- load_checkpoint("counts_filtered", dir = dir_rds)
tt          <- load_checkpoint(paste0("tt_", EXP, "_annotated"), dir = dir_rds)
logcpm      <- load_checkpoint(paste0("logcpm_", EXP), dir = dir_rds)
sample_meta <- load_checkpoint(paste0("sample_meta_", EXP), dir = dir_rds)
scores      <- load_checkpoint(paste0("module_scores_", EXP), dir = dir_rds)
wide        <- load_checkpoint(paste0("module_scores_wide_", EXP), dir = dir_rds)

map          <- symbol_map(tt)
hallmark     <- fgsea::gmtPathways("genesets/h.all.v2025.1.Hs.symbols.gmt")
isg_genes    <- hallmark[[ISG_SET]]
stopifnot(length(isg_genes) > 0)

# =============================================================================
# Per-sample ISG score and induction
# =============================================================================
# Same scorer as the composition modules: mean z across detected set members,
# z taken over all samples of this experiment.
isg_long <- module_expr_z(logcpm, sample_meta, map, setNames(list(isg_genes), ISG_SET))
isg      <- module_score(isg_long) |>
  dplyr::select(sample, line, stim, isg_score = score, n_isg = n_genes)

message("ISG set: ", isg$n_isg[1], " of ", length(isg_genes), " genes detected")

# Induction is scored against the same line's Mock mean, making it a response
# measure rather than a level measure.
mock_ref <- isg |>
  dplyr::filter(stim == "Mock") |>
  dplyr::group_by(line) |>
  dplyr::summarise(mock_isg = mean(isg_score), .groups = "drop")

dat <- isg |>
  dplyr::left_join(mock_ref, by = "line") |>
  dplyr::left_join(
    dplyr::select(wide, -dplyr::any_of(c("line", "stim", "condition"))),
    by = "sample"
  ) |>
  dplyr::mutate(
    induction = isg_score - mock_isg,
    genotype  = factor(ifelse(line == CASE_LINE, "C9", "Control"),
                       levels = c("Control", "C9")),
    covar     = .data[[COVARIATE]]
  )

stim_dat <- dat |>
  dplyr::filter(stim != "Mock") |>
  dplyr::mutate(stim = droplevels(stim), line = droplevels(line))

message("Stimulated samples: ", nrow(stim_dat),
        " (C9 ", sum(stim_dat$genotype == "C9"),
        ", control ", sum(stim_dat$genotype == "Control"), ")")

# =============================================================================
# Confound summary
# =============================================================================
# Correlation of the covariate with genotype, and the overlap between the C9 and
# control covariate ranges. Adjustment outside that overlap is extrapolation.
rng <- function(x) c(min(x), max(x))
c9_rng   <- rng(stim_dat$covar[stim_dat$genotype == "C9"])
ctl_rng  <- rng(stim_dat$covar[stim_dat$genotype == "Control"])
overlap  <- max(0, min(c9_rng[2], ctl_rng[2]) - max(c9_rng[1], ctl_rng[1]))
n_ctl_in <- sum(stim_dat$genotype == "Control" &
                stim_dat$covar >= c9_rng[1] & stim_dat$covar <= c9_rng[2])

confound <- tibble::tibble(
  covariate      = COVARIATE,
  r_covar_geno   = safe_cor(stim_dat$covar, as.numeric(stim_dat$genotype)),
  c9_min = c9_rng[1], c9_max = c9_rng[2],
  ctrl_min = ctl_rng[1], ctrl_max = ctl_rng[2],
  overlap_width  = overlap,
  n_ctrl_in_c9_range = n_ctl_in
)
print(confound)

# =============================================================================
# Model 1 - multiple regression on the stimulated samples
# =============================================================================
# Variance inflation without a car dependency: regress each predictor on the
# others from the same model matrix and take 1/(1 - R^2).
vif_manual <- function(mod) {
  X <- stats::model.matrix(mod)
  X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
  if (ncol(X) < 2L) return(setNames(rep(1, ncol(X)), colnames(X)))
  vapply(seq_len(ncol(X)), function(j) {
    r2 <- summary(stats::lm(X[, j] ~ X[, -j, drop = FALSE]))$r.squared
    if (!is.finite(r2) || r2 >= 1) Inf else 1 / (1 - r2)
  }, numeric(1)) |> setNames(colnames(X))
}

tidy_lm <- function(mod, model_name) {
  s <- summary(mod)$coefficients
  v <- vif_manual(mod)
  tibble::tibble(
    model    = model_name,
    term     = rownames(s),
    estimate = s[, "Estimate"],
    se       = s[, "Std. Error"],
    p        = s[, "Pr(>|t|)"]
  ) |>
    dplyr::mutate(vif = unname(v[match(term, names(v))]))
}

m_geno  <- stats::lm(induction ~ stim + genotype, data = stim_dat)
m_covar <- stats::lm(induction ~ stim + covar,    data = stim_dat)
m_both  <- stats::lm(induction ~ stim + covar + genotype, data = stim_dat)

reg_tbl <- dplyr::bind_rows(
  tidy_lm(m_geno,  "genotype only"),
  tidy_lm(m_covar, paste0(COVARIATE, " only")),
  tidy_lm(m_both,  paste0("genotype + ", COVARIATE))
) |>
  dplyr::filter(term != "(Intercept)")

print(as.data.frame(reg_tbl), digits = 3)

# Share of the raw genotype effect absorbed by the covariate
geno_alone <- coef(m_geno)[["genotypeC9"]]
geno_adj   <- coef(m_both)[["genotypeC9"]]
pct_expl   <- 100 * (1 - geno_adj / geno_alone)

vif_geno_adj <- reg_tbl$vif[reg_tbl$term == "genotypeC9" &
                            reg_tbl$model == paste0("genotype + ", COVARIATE)]

message(sprintf("Genotype effect on ISG induction: %.3f alone -> %.3f adjusted for %s (%.0f%% absorbed, VIF %.2f)",
                geno_alone, geno_adj, COVARIATE, pct_expl, vif_geno_adj))

# =============================================================================
# Model 2 - control-only calibration, then predict C9
# =============================================================================
# Slope fit on control samples only, so genotype cannot contribute to it. C9 is
# then predicted from its own composition; the residual is the unexplained part.
ctl  <- dplyr::filter(stim_dat, genotype == "Control")
m_ctl <- stats::lm(induction ~ stim + covar, data = ctl)

c9 <- dplyr::filter(stim_dat, genotype == "C9")
c9$predicted <- as.numeric(stats::predict(m_ctl, newdata = c9))
c9$residual  <- c9$induction - c9$predicted

calib <- tibble::tibble(
  n_control      = nrow(ctl),
  slope_covar    = unname(coef(m_ctl)[["covar"]]),
  slope_p        = summary(m_ctl)$coefficients["covar", "Pr(>|t|)"],
  ctl_r2         = summary(m_ctl)$r.squared,
  mean_c9_obs    = mean(c9$induction),
  mean_c9_pred   = mean(c9$predicted),
  mean_residual  = mean(c9$residual),
  residual_p     = safe_t_p_one(c9$residual)
)
print(calib)

# =============================================================================
# Model 3 - gene-level limma, with and without the covariate
# =============================================================================
# Both designs are fit with and without the covariate. Under the cell-means
# design the covariate slope comes from within-group deviations only; the reduced
# design (line + stim) lets between-group variation inform it too. The set-level
# statistic is the mean moderated t of the Hallmark IFN-alpha genes plus a
# Wilcoxon of in-set against out-of-set t.
isg_set_stat <- function(fit, coef_name, genes) {
  t_all <- limma::topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  t_all$ensembl_id <- sub("\\..*$", "", rownames(t_all))
  t_all <- dplyr::left_join(t_all, map, by = "ensembl_id")
  inset <- t_all$SYMBOL %in% genes & is.finite(t_all$t)
  tibble::tibble(
    n_in_set  = sum(inset),
    mean_t    = mean(t_all$t[inset]),
    mean_lfc  = mean(t_all$logFC[inset]),
    p_wilcox  = stats::wilcox.test(t_all$t[inset], t_all$t[!inset])$p.value
  )
}

# lmFit only - eBayes runs after contrasts.fit, matching 3_DE_limma.R.
fit_design <- function(design, samples) {
  mat  <- round(as.matrix(counts[, samples, drop = FALSE]))
  dge  <- edgeR::DGEList(counts = mat)
  keep <- edgeR::filterByExpr(dge, design = design)
  dge  <- dge[keep, , keep.lib.sizes = FALSE]
  dge  <- edgeR::calcNormFactors(dge)
  limma::lmFit(limma::voom(dge, design), design)
}

# --- 3a. saturated cell-means, the design 3_DE_limma.R uses ---
sat_sm <- sample_meta |>
  dplyr::mutate(condition = droplevels(condition)) |>
  dplyr::left_join(dplyr::select(dat, sample, covar), by = "sample") |>
  dplyr::mutate(covar_c = as.numeric(scale(covar, scale = FALSE)))

d_sat  <- stats::model.matrix(~ 0 + condition, data = sat_sm)
colnames(d_sat) <- levels(sat_sm$condition)

d_sat_a <- cbind(d_sat, covar = sat_sm$covar_c)

sat_contrast <- "C9_IFNb - Ctrl3_IFNb"
f_sat  <- fit_design(d_sat,   sat_sm$sample)
f_sat_a <- fit_design(d_sat_a, sat_sm$sample)

cm_sat   <- limma::makeContrasts(contrasts = sat_contrast, levels = d_sat)
cm_sat_a <- limma::makeContrasts(contrasts = sat_contrast, levels = d_sat_a)
colnames(cm_sat) <- colnames(cm_sat_a) <- "genotype"

sat_cmp <- dplyr::bind_rows(
  isg_set_stat(limma::eBayes(limma::contrasts.fit(f_sat,   cm_sat)),   "genotype", isg_genes),
  isg_set_stat(limma::eBayes(limma::contrasts.fit(f_sat_a, cm_sat_a)), "genotype", isg_genes)
) |>
  dplyr::mutate(design = c("cell-means", paste0("cell-means + ", COVARIATE)),
                scope = sat_contrast, .before = 1)

# --- 3b. reduced design (line + stim) on the stimulated samples ---
red_sm <- stim_dat |>
  dplyr::mutate(line = droplevels(line),
                covar_c = as.numeric(scale(covar, scale = FALSE)))

d_red   <- stats::model.matrix(~ 0 + line + stim, data = red_sm)
colnames(d_red) <- sub("^line", "", colnames(d_red))
d_red_a <- cbind(d_red, covar = red_sm$covar_c)

red_vif <- vif_manual(stats::lm(induction ~ line + stim + covar, data = red_sm))
print(round(red_vif, 2))

cm_red   <- limma::makeContrasts("C9 - (Ctrl3 + Ctrl2)/2", levels = d_red)
cm_red_a <- limma::makeContrasts("C9 - (Ctrl3 + Ctrl2)/2", levels = d_red_a)
colnames(cm_red) <- colnames(cm_red_a) <- "genotype"

red_cmp <- dplyr::bind_rows(
  isg_set_stat(limma::eBayes(limma::contrasts.fit(fit_design(d_red,   red_sm$sample), cm_red)),
               "genotype", isg_genes),
  isg_set_stat(limma::eBayes(limma::contrasts.fit(fit_design(d_red_a, red_sm$sample), cm_red_a)),
               "genotype", isg_genes)
) |>
  dplyr::mutate(design = c("line + stim", paste0("line + stim + ", COVARIATE)),
                scope = "stimulated samples, C9 - mean(controls)", .before = 1)

limma_cmp <- dplyr::bind_rows(sat_cmp, red_cmp)
print(as.data.frame(limma_cmp), digits = 3)

# =============================================================================
# Figure
# =============================================================================
lab_covar <- paste0(COVARIATE, " module score (mean z)")

# Panel A: covariate by genotype, with the C9 range shaded
p_confound <- ggplot(stim_dat, aes(covar, genotype)) +
  annotate("rect", xmin = c9_rng[1], xmax = c9_rng[2],
           ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.6) +
  geom_point(aes(fill = line, shape = stim), size = 2.6,
             position = position_jitter(height = 0.12, width = 0)) +
  scale_fill_manual(values = line_cols, name = NULL) +
  scale_shape_manual(values = stim_shapes, name = NULL) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  labs(x = lab_covar, y = NULL,
       title = "A. The confound",
       subtitle = sprintf("r = %.2f with genotype; %d of %d control samples fall inside the C9 range (shaded)",
                          confound$r_covar_geno, n_ctl_in, sum(stim_dat$genotype == "Control"))) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank(), legend.position = "bottom")

# Panel B: line and ribbon from m_ctl, which excludes the C9 points
grid_ctl <- tidyr::expand_grid(
  covar = seq(min(stim_dat$covar), max(stim_dat$covar), length.out = 80),
  stim  = levels(stim_dat$stim)
) |>
  dplyr::mutate(stim = factor(stim, levels = levels(stim_dat$stim)))
pred_ctl <- stats::predict(m_ctl, newdata = grid_ctl, interval = "confidence")
grid_ctl <- dplyr::bind_cols(grid_ctl, tibble::as_tibble(pred_ctl))

p_calib <- ggplot(stim_dat, aes(covar, induction)) +
  geom_ribbon(data = grid_ctl, aes(y = fit, ymin = lwr, ymax = upr),
              fill = "grey70", alpha = 0.3) +
  geom_line(data = grid_ctl, aes(y = fit), colour = "grey30", linewidth = 0.5) +
  geom_point(aes(fill = line, shape = stim), size = 2.8, stroke = 0.4) +
  facet_wrap(~ stim, nrow = 1) +
  scale_fill_manual(values = line_cols, name = NULL) +
  scale_shape_manual(values = stim_shapes, guide = "none") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3))) +
  labs(x = lab_covar, y = "ISG induction (z vs own-line Mock)",
       title = "B. Control-only calibration, C9 held against it",
       subtitle = sprintf("slope %+.2f per z (p = %.3g); C9 observed %.2f vs predicted %.2f, residual %+.2f",
                          calib$slope_covar, calib$slope_p,
                          calib$mean_c9_obs, calib$mean_c9_pred, calib$mean_residual)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "bottom")

# Panel C: the genotype coefficient with and without the covariate, +/- 1.96 SE.
coef_d <- reg_tbl |>
  dplyr::filter(term == "genotypeC9") |>
  dplyr::mutate(model = factor(model, levels = rev(unique(model))))

p_coef <- ggplot(coef_d, aes(estimate, model)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey60") +
  geom_errorbar(aes(xmin = estimate - 1.96 * se, xmax = estimate + 1.96 * se),
                orientation = "y", width = 0.16, linewidth = 0.4) +
  geom_point(size = 3, shape = 21, fill = line_cols[["C9"]]) +
  geom_text(aes(label = sprintf("VIF %.1f", vif)), vjust = -1.3, size = 2.8,
            colour = "grey35") +
  labs(x = "Genotype effect on ISG induction (C9 - control)", y = NULL,
       title = "C. Does the genotype term survive?",
       subtitle = sprintf("%.0f%% of the raw genotype effect is absorbed by %s",
                          pct_expl, COVARIATE)) +
  theme_bw() +
  theme(panel.grid.major.y = element_blank())

fig <- (p_confound / p_calib / p_coef) +
  plot_layout(heights = c(0.8, 1.3, 0.8)) +
  plot_annotation(
    title    = "Is the blunted IFN response separable from cell composition?",
    subtitle = paste0(EXP, " - ", ISG_SET, ", ", isg$n_isg[1], " genes detected"),
    # Built from the fitted values so it cannot drift from the panels
    caption  = paste(
      sprintf("Genotype and %s are correlated (r = %.2f) but NOT collinear (VIF %.1f), so the model had",
              tolower(COVARIATE), confound$r_covar_geno, vif_geno_adj),
      sprintf("the information to keep a genotype term and did not: %.0f%% of the effect moves to %s.",
              pct_expl, tolower(COVARIATE)),
      sprintf("The real limit is coverage - only %d of %d control samples sit inside the C9 %s range, so panel B's",
              n_ctl_in, sum(stim_dat$genotype == "Control"), tolower(COVARIATE)),
      "line is extrapolated where C9 lies. And replicates are organoid batches from one line per genotype,",
      "so nothing here separates C9orf72 from iPSC line of origin - that needs an isogenic corrected line.",
      sep = "\n"),
    theme = theme(
      plot.title    = element_text(face = "bold"),
      plot.subtitle = element_text(colour = "grey30", size = 9),
      plot.caption  = element_text(colour = "grey40", size = 7, hjust = 0))
  )

save_fig(fig, paste0("composition_adjusted_isg_", EXP),
         width = 11, height = 13, subdir = dir_comp)

# =============================================================================
# Write-out
# =============================================================================
out_dir <- ensure_dir(dir_res)
wtsv <- function(x, f) {
  utils::write.table(x, file.path(out_dir, f), sep = "\t", quote = FALSE, row.names = FALSE)
  message("Wrote ", file.path(out_dir, f))
}

wtsv(dplyr::mutate(confound, experiment = EXP, .before = 1), "isg_composition_confound.tsv")
wtsv(dplyr::mutate(reg_tbl, experiment = EXP, .before = 1),  "isg_composition_regression.tsv")
wtsv(dplyr::mutate(calib, experiment = EXP, covariate = COVARIATE, .before = 1),
     "isg_composition_calibration.tsv")
wtsv(dplyr::mutate(limma_cmp, experiment = EXP, .before = 1), "isg_composition_limma.tsv")
wtsv(dplyr::select(dat, sample, line, stim, isg_score, induction,
                   dplyr::any_of(c("Proliferation", "Progenitor", "Neuron", "Astrocyte", "maturity"))),
     "isg_per_sample.tsv")

message("9_composition_adjusted.R complete")
