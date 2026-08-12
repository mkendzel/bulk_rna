# =============================================
# 2 - limma-voom differential expression
# =============================================
# Loads `counts_filtered` from 1_mavs_QC.R and runs limma-voom DE for
# three comparisons:
#   - Spleen D0 vs D7 (WT and MAVSKO)
#   - Brain  D0 vs D7 (WT and MAVSKO)
#   - Brain vs Spleen at D7 (Spleen = reference)
# Saves voom objects and (symbol-annotated) topTable result lists.
# Output checkpoints feed 3_mavs_GSEA.R, 4_mavs_TF_activity.R, and 5_mavs_KEGG.R.

# ---- Libraries ----
library(tidyverse)
library(edgeR)
library(limma)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

# =============================================
# --------------- Spleen D0 vs D7 --------------
# =============================================
# subset counts to spleen samples only
counts <- load_checkpoint(
   "counts_filtered",
   dir = "projects/mavs_mito/data/r_objects"
)


sample_names <- colnames(counts)
colData <- data.frame(sample = sample_names) |>
  dplyr::mutate(
    condition = sub("^[0-9]+_", "", sample)
  )
rownames(colData) <- colData$sample

spleen_mock_samples <- sample_names[
   (grepl("Spleen", sample_names) | grepl("_d0$", sample_names)) &
     !grepl("Brain", sample_names)
]

counts_spleen <- counts[, spleen_mock_samples]


spleen_d0_d7_samples <- colnames(counts_spleen)[
  grepl("_d0$", colnames(counts_spleen)) | grepl("_D7$", colnames(counts_spleen))
]

counts_spleen_d0_d7 <- counts_spleen[, spleen_d0_d7_samples]

genotype <- ifelse(grepl("MAVSko|MAVSKO", spleen_d0_d7_samples), "MAVSKO", "WT")
day <- ifelse(grepl("_d0$", spleen_d0_d7_samples), "D0", "D7")

colData_spleen <- data.frame(
   sample = spleen_d0_d7_samples,
   genotype = genotype,
   day = day,
   group = factor(paste(genotype, day, sep = "_"))
)

save_checkpoint(colData_spleen, "colData_spleen", dir = "projects/mavs_mito/data/r_objects",
                lines = c(19:60),
                notes = "Meta data for spleen D0 vs D7 samples")

# Limma-voom DE: WT and MAVSKO, naive (D0) vs day 7 post-infection

mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")

# Check rowname format before stripping versions
head(rownames(counts_spleen_d0_d7))

# counts_spleen_d0_d7 and colData_spleen built previously
stopifnot(identical(colnames(counts_spleen_d0_d7), colData_spleen$sample))

colData_spleen$group <- factor(colData_spleen$group,
  levels = c("WT_D0", "MAVSKO_D0", "WT_D7", "MAVSKO_D7")
)
table(colData_spleen$group)

# Create DGEList and filter
dge <- DGEList(counts = as.matrix(counts_spleen_d0_d7), group = colData_spleen$group)
keep <- filterByExpr(dge, group = colData_spleen$group) # keep genes expressed in at least half of the samples in any group
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cat("Genes before filtering:", nrow(counts_spleen_d0_d7), "\n")
cat("Genes after filtering:", nrow(dge), "\n")

# Design matrix and voom
design <- model.matrix(~ 0 + group, data = colData_spleen)
colnames(design) <- levels(colData_spleen$group)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# Contrasts: within-genotype day effects, direct genotype comparison, interaction
contr_matrix <- makeContrasts(
  WT_D7_vs_D0 = WT_D7 - WT_D0,
  MAVSKO_D7_vs_D0 = MAVSKO_D7 - MAVSKO_D0,
  WT_vs_MAVSKO_D7 = WT_D7 - MAVSKO_D7,
  interaction = (WT_D7 - WT_D0) - (MAVSKO_D7 - MAVSKO_D0),
  levels = design
)

fit2 <- contrasts.fit(fit, contr_matrix)
fit2 <- eBayes(fit2)

res_voom_spleen_d0_d7 <- setNames(
  lapply(colnames(contr_matrix), function(coef) {
    topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  }),
  colnames(contr_matrix)
)

# Quick sanity check on each contrast
lapply(res_voom_spleen_d0_d7, function(df) sum(df$adj.P.Val < 0.05))

# Annotate with gene symbols
res_voom_spleen_d0_d7_annotated <- lapply(res_voom_spleen_d0_d7, function(df) {
  df %>%
    rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) %>%
    left_join(
      mouse_gene_map %>% distinct(ENSEMBL, .keep_all = TRUE),
      by = c("ensembl_id" = "ENSEMBL")
    ) %>%
    as_tibble()
})

# Check annotation join coverage
lapply(res_voom_spleen_d0_d7_annotated, function(df) sum(is.na(df$SYMBOL)))

# Save
save_checkpoint(v, "voom_spleen_d0_d7",
  notes = "voom object for spleen D0 vs D7, WT and MAVSKO. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(21:135)
)

save_checkpoint(res_voom_spleen_d0_d7, "res_voom_spleen_d0_d7",
  notes = "limma-voom DE results for spleen D0 vs D7; named list of topTables (within-genotype day effects, direct genotype comparison at D7, interaction). Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(21:135)
)

save_checkpoint(res_voom_spleen_d0_d7_annotated, "res_voom_spleen_d0_d7_annotated",
  notes = "Gene Symbol annotated limma-voom DE results for spleen D0 vs D7. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(21:135)
)

# =============================================
# --------------- Brain D0 vs D7 ---------------
# =============================================
counts <- load_checkpoint(
  "counts_filtered",
  dir = "projects/mavs_mito/data/r_objects"
)

sample_names <- colnames(counts)

brain_mock_samples <- sample_names[
  (grepl("Brain", sample_names) | grepl("_d0$", sample_names)) &
    !grepl("Spleen", sample_names)
]

counts_brain <- counts[, brain_mock_samples]

brain_d0_d7_samples <- colnames(counts_brain)[
  grepl("_d0$", colnames(counts_brain)) | grepl("_D7$", colnames(counts_brain))
]

counts_brain_d0_d7 <- counts_brain[, brain_d0_d7_samples]

genotype <- ifelse(grepl("MAVSko|MAVSKO", brain_d0_d7_samples), "MAVSKO", "WT")
day <- ifelse(grepl("_d0$", brain_d0_d7_samples), "D0", "D7")

colData_brain <- data.frame(
  sample = brain_d0_d7_samples,
  genotype = genotype,
  day = day,
  group = factor(paste(genotype, day, sep = "_"))
)

save_checkpoint(colData_brain, "colData_brain", dir = "projects/mavs_mito/data/r_objects",
                lines = c(163:196),
                notes = "Meta data for brain D0 vs D7 samples")

# Limma-voom DE: WT and MAVSKO, naive (D0) vs day 7 post-infection

mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")

head(rownames(counts_brain_d0_d7))

stopifnot(identical(colnames(counts_brain_d0_d7), colData_brain$sample))

colData_brain$group <- factor(colData_brain$group,
  levels = c("WT_D0", "MAVSKO_D0", "WT_D7", "MAVSKO_D7")
)
table(colData_brain$group)

# Create DGEList and filter
dge <- DGEList(counts = as.matrix(counts_brain_d0_d7), group = colData_brain$group)
keep <- filterByExpr(dge, group = colData_brain$group)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

cat("Genes before filtering:", nrow(counts_brain_d0_d7), "\n")
cat("Genes after filtering:", nrow(dge), "\n")

# Design matrix and voom
design <- model.matrix(~ 0 + group, data = colData_brain)
colnames(design) <- levels(colData_brain$group)
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# Contrasts: within-genotype day effects, direct genotype comparison, interaction
contr_matrix <- makeContrasts(
  WT_D7_vs_D0 = WT_D7 - WT_D0,
  MAVSKO_D7_vs_D0 = MAVSKO_D7 - MAVSKO_D0,
  WT_vs_MAVSKO_D7 = WT_D7 - MAVSKO_D7,
  interaction = (WT_D7 - WT_D0) - (MAVSKO_D7 - MAVSKO_D0),
  levels = design
)

fit2 <- contrasts.fit(fit, contr_matrix)
fit2 <- eBayes(fit2)

res_voom_brain_d0_d7 <- setNames(
  lapply(colnames(contr_matrix), function(coef) {
    topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  }),
  colnames(contr_matrix)
)

# Quick sanity check on each contrast
lapply(res_voom_brain_d0_d7, function(df) sum(df$adj.P.Val < 0.05))

# Annotate with gene symbols
res_voom_brain_d0_d7_annotated <- lapply(res_voom_brain_d0_d7, function(df) {
  df %>%
    rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) %>%
    left_join(
      mouse_gene_map %>% distinct(ENSEMBL, .keep_all = TRUE),
      by = c("ensembl_id" = "ENSEMBL")
    ) %>%
    as_tibble()
})

# Check annotation join coverage
lapply(res_voom_brain_d0_d7_annotated, function(df) sum(is.na(df$SYMBOL)))

# Save
save_checkpoint(v, "voom_brain_d0_d7",
  notes = "voom object for brain D0 vs D7, WT and MAVSKO. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(163:249)
)

save_checkpoint(res_voom_brain_d0_d7, "res_voom_brain_d0_d7",
  notes = "limma-voom DE results for brain D0 vs D7; named list of topTables (within-genotype day effects, direct genotype comparison at D7, interaction). Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(163:249)
)

save_checkpoint(res_voom_brain_d0_d7_annotated, "res_voom_brain_d0_d7_annotated",
  notes = "Gene Symbol annotated limma-voom DE results for brain D0 vs D7. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(163:249)
)

# =============================================
# ------ Brain D7 vs Spleen D7 (ref = Spleen) ------
# =============================================
counts <- load_checkpoint("counts_filtered", dir = "projects/mavs_mito/data/r_objects")
sample_names <- colnames(counts)

# All D7 spleen + brain samples, both genotypes
d7_tissue_samples <- sample_names[
  grepl("_D7$", sample_names) &
    (grepl("Spleen", sample_names) | grepl("Brain", sample_names))
]
counts_d7_tissue <- counts[, d7_tissue_samples]

genotype <- ifelse(grepl("MAVSko|MAVSKO", d7_tissue_samples), "MAVSKO", "WT")
tissue   <- ifelse(grepl("Brain", d7_tissue_samples), "Brain", "Spleen")

colData_d7_tissue <- data.frame(
  sample   = d7_tissue_samples,
  genotype = genotype,
  tissue   = tissue,
  # Spleen listed first so it is the reference level within each genotype
  group    = factor(paste(genotype, tissue, sep = "_"),
                    levels = c("WT_Spleen", "WT_Brain", "MAVSKO_Spleen", "MAVSKO_Brain"))
)
stopifnot(identical(colnames(counts_d7_tissue), colData_d7_tissue$sample))
table(colData_d7_tissue$group)   # expect 5 each

save_checkpoint(colData_d7_tissue, "colData_d7_tissue", dir = "projects/mavs_mito/data/r_objects",
                lines = c(268:290),
                notes = "Meta data for Brain vs Spleen D7 tissue comparison (both genotypes)")

# Limma-voom DE: Brain vs Spleen at day 7, WT and MAVSKO separately

mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")

head(rownames(counts_d7_tissue))

# DGEList + filter (group-aware), voom, fit
dge  <- DGEList(counts = as.matrix(counts_d7_tissue), group = colData_d7_tissue$group)
keep <- filterByExpr(dge, group = colData_d7_tissue$group)
dge  <- calcNormFactors(dge[keep, , keep.lib.sizes = FALSE])

cat("Genes before filtering:", nrow(counts_d7_tissue), "\n")
cat("Genes after filtering:", nrow(dge), "\n")

design <- model.matrix(~ 0 + group, data = colData_d7_tissue)
colnames(design) <- levels(colData_d7_tissue$group)
v   <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)

# Brain - Spleen within each genotype (Spleen = reference)
contr_matrix <- makeContrasts(
  WT_Brain_vs_Spleen     = WT_Brain     - WT_Spleen,
  MAVSKO_Brain_vs_Spleen = MAVSKO_Brain - MAVSKO_Spleen,
  levels = design
)
fit2 <- eBayes(contrasts.fit(fit, contr_matrix))

res_voom_d7_tissue <- setNames(
  lapply(colnames(contr_matrix), function(coef) {
    topTable(fit2, coef = coef, number = Inf, sort.by = "none")
  }),
  colnames(contr_matrix)
)

# Quick sanity check on each contrast
lapply(res_voom_d7_tissue, function(df) sum(df$adj.P.Val < 0.05))

# Annotate with gene symbols (same join used elsewhere)
res_voom_d7_tissue_annotated <- lapply(res_voom_d7_tissue, function(df) {
  df %>%
    rownames_to_column("ensembl_id") %>%
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) %>%
    left_join(
      mouse_gene_map %>% distinct(ENSEMBL, .keep_all = TRUE),
      by = c("ensembl_id" = "ENSEMBL")
    ) %>%
    as_tibble()
})

# Check annotation join coverage
lapply(res_voom_d7_tissue_annotated, function(df) sum(is.na(df$SYMBOL)))

# Save
save_checkpoint(v, "voom_d7_tissue",
  notes = "voom object: Brain vs Spleen at D7, 4-level genotype_tissue model. Default filtering, min count = 10.",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(268:326)
)

save_checkpoint(res_voom_d7_tissue, "res_voom_d7_tissue",
  notes = "limma-voom DE: Brain - Spleen at D7 within WT and within MAVSKO (Spleen = reference). Named list of topTables.",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(268:326)
)

save_checkpoint(res_voom_d7_tissue_annotated, "res_voom_d7_tissue_annotated",
  notes = "Gene Symbol annotated Brain - Spleen D7 DE (WT and MAVSKO).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(268:326)
)
