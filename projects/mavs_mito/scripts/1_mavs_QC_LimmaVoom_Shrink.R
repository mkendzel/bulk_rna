# ---- Libraries ----
library(stringr)
library(tidyverse)
library(DESeq2)
library(ggvenn)
library(grid)
library(readxl)
library(edgeR)
library(limma)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

# ----- Import data ----
## Primary import: expression matrix (.txt)
expr <- read.delim(
  "projects/mavs_mito/data/raw/All samples.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(expr) <- expr$gene_id
expr <- expr[, -1]

## Alternative import: subset of expression matrix (.txt)
sub_expr <- read.delim(
  "projects/mavs_mito/data/raw/527_EU_Suthar_60RNAseq_4_RawCounts_wGeneSymbols.txt",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

#remove first four rows
sub_expr <- sub_expr[-(1:4), ]

#Rename rownames to gene_id
rownames(sub_expr) <- sub_expr$gene_id
sub_expr <- sub_expr[, -1]


# ---- Explore samples ----
## 1) Count unique samples in each data set
# "all data set"
get_suffix <- function(cols) sub("^[0-9]+_", "", cols)

expr_suffixes <- get_suffix(colnames(expr))

suffix_counts <- as.data.frame(table(expr_suffixes), stringsAsFactors = FALSE)
names(suffix_counts) <- c("suffix", "n")
suffix_counts <- suffix_counts[order(-suffix_counts$n), ]
rownames(suffix_counts) <- NULL

suffix_counts

# "subset data set"
get_suffix_drop2 <- function(cols) sub("^[^_]*_[^_]*_", "", cols)

sub_expr_suffixes <- get_suffix_drop2(colnames(sub_expr))

sub_suffix_counts <- as.data.frame(table(sub_expr_suffixes), stringsAsFactors = FALSE)
names(sub_suffix_counts) <- c("suffix", "n")
sub_suffix_counts <- sub_suffix_counts[order(-sub_suffix_counts$n), ]
rownames(sub_suffix_counts) <- NULL

sub_suffix_counts

# Table of suffixes in both data sets
combined_counts <- merge(
  suffix_counts, sub_suffix_counts,
  by = "suffix", all = TRUE,
  suffixes = c("_expr", "_sub")
)
names(combined_counts) <- c("suffix", "n_expr", "n_sub")
combined_counts[is.na(combined_counts)] <- 0
combined_counts <- combined_counts[order(-combined_counts$n_expr, -combined_counts$n_sub), ]
rownames(combined_counts) <- NULL

knitr::kable(combined_counts, format = "markdown")

# Notes:
# Missing 3 Spleen Day 30 samples in the "All data set" stored in the object "expr"
# Solution: Include the missing samples to "All data" from the "subset data set"
# Pass these through QC to see why they were removed

# ---- Compare genes and gene counts between data sets ----
# Extract leading sample number from column names
get_sample_num <- function(cols) as.integer(sub("^([0-9]+).*", "\\1", cols))

expr_num <- get_sample_num(colnames(expr))
sub_num <- get_sample_num(colnames(sub_expr))

expr_lookup <- data.frame(col = colnames(expr), num = expr_num)
expr_lookup <- expr_lookup[!is.na(expr_lookup$num), ]

sub_lookup <- data.frame(col = colnames(sub_expr), num = sub_num)
sub_lookup <- sub_lookup[!is.na(sub_lookup$num), ]

## 1) Do the two data sets share the same samples? ----
samples_only_in_expr <- sort(setdiff(expr_lookup$num, sub_lookup$num))
samples_only_in_sub <- sort(setdiff(sub_lookup$num, expr_lookup$num))
samples_dup_in_expr <- sort(unique(expr_lookup$num[duplicated(expr_lookup$num)]))
samples_dup_in_sub <- sort(unique(sub_lookup$num[duplicated(sub_lookup$num)]))

same_samples <- length(samples_only_in_expr) == 0 &&
  length(samples_only_in_sub) == 0 &&
  length(samples_dup_in_expr) == 0 &&
  length(samples_dup_in_sub) == 0

sample_check <- list(
  same_samples = same_samples,
  samples_only_in_expr = samples_only_in_expr,
  samples_only_in_sub = samples_only_in_sub,
  duplicated_sample_numbers_in_expr = samples_dup_in_expr,
  duplicated_sample_numbers_in_sub = samples_dup_in_sub
)
sample_check

# Sample numbers safe to compare: present in both, not duplicated in either
shared_nums <- sort(setdiff(
  intersect(expr_lookup$num, sub_lookup$num),
  union(samples_dup_in_expr, samples_dup_in_sub)
))

## 2) Within shared samples, how many genes are shared? ----
shared_genes <- intersect(rownames(expr), rownames(sub_expr))

gene_check <- list(
  n_genes_expr = nrow(expr),
  n_genes_sub = nrow(sub_expr),
  n_genes_shared = length(shared_genes),
  genes_only_in_expr = setdiff(rownames(expr), rownames(sub_expr)),
  genes_only_in_sub = setdiff(rownames(sub_expr), rownames(expr))
)
gene_check

## 3) Do shared genes have identical counts? ----
# not present as comparable conditions in sub_expr
cols_to_remove <- c(
  "1_WT_d0", "2_WT_d0", "3_WT_d0", "4_WT_d0",
  "5_MAVSKO_d0", "6_MAVSKO_d0", "7_MAVSKO_d0", "8_MAVSKO_d0"
)

expr_filtered <- expr[, !(colnames(expr) %in% cols_to_remove)]

# Rebuild sample lookup tables using the filtered expr
expr_num <- get_sample_num(colnames(expr_filtered))
expr_lookup <- data.frame(col = colnames(expr_filtered), num = expr_num)
expr_lookup <- expr_lookup[!is.na(expr_lookup$num), ]

sub_num <- get_sample_num(colnames(sub_expr))
sub_lookup <- data.frame(col = colnames(sub_expr), num = sub_num)
sub_lookup <- sub_lookup[!is.na(sub_lookup$num), ]

# Recheck for duplicates/unmatched samples after filtering
samples_only_in_expr <- sort(setdiff(expr_lookup$num, sub_lookup$num))
samples_only_in_sub <- sort(setdiff(sub_lookup$num, expr_lookup$num))
samples_dup_in_expr <- sort(unique(expr_lookup$num[duplicated(expr_lookup$num)]))
samples_dup_in_sub <- sort(unique(sub_lookup$num[duplicated(sub_lookup$num)]))

shared_nums <- sort(setdiff(
  intersect(expr_lookup$num, sub_lookup$num),
  union(samples_dup_in_expr, samples_dup_in_sub)
))

# Shared genes between filtered expr and sub_expr
shared_genes <- intersect(rownames(expr_filtered), rownames(sub_expr))

compare_one_sample <- function(n) {
  expr_col <- expr_lookup$col[expr_lookup$num == n]
  sub_col <- sub_lookup$col[sub_lookup$num == n]

  x <- as.numeric(expr_filtered[shared_genes, expr_col])
  y <- as.numeric(sub_expr[shared_genes, sub_col])
  is_match <- x == y

  data.frame(
    number = n,
    expr_col = expr_col,
    sub_col = sub_col,
    n_genes = length(shared_genes),
    n_mismatch = sum(!is_match, na.rm = TRUE),
    pct_mismatch = round(mean(!is_match, na.rm = TRUE) * 100, 2)
  )
}

sample_comparison <- do.call(rbind, lapply(shared_nums, compare_one_sample))
sample_comparison <- sample_comparison[order(sample_comparison$number), ]
rownames(sample_comparison) <- NULL

sample_comparison
mismatched_samples <- sample_comparison[sample_comparison$n_mismatch > 0, ]
mismatched_samples

# Notes: There are no missmatched counts. The data sets are identical except for the missing Spleen day 30s

# ---- Edit expr to include missing samples ----
# Remove all Spleen_D30 columns from expr
expr_no_d30 <- expr[, !grepl("Spleen_D30", colnames(expr))]

# Add all Spleen_D30 columns from sub_expr
d30_cols_sub <- colnames(sub_expr)[grepl("Spleen_D30", colnames(sub_expr))]
expr_updated <- cbind(expr_no_d30, sub_expr[rownames(expr_no_d30), d30_cols_sub])

#Check
get_suffix <- function(cols) sub("^[0-9]+_", "", cols)

expr_suffixes <- get_suffix(colnames(expr_updated))

suffix_counts <- as.data.frame(table(expr_suffixes), stringsAsFactors = FALSE)
names(suffix_counts) <- c("suffix", "n")
suffix_counts <- suffix_counts[order(-suffix_counts$n), ]
rownames(suffix_counts) <- NULL

suffix_counts 

# Save checkpoint
save_checkpoint(
  expr_updated,
  "all_samples_updated",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(12:99),
  notes = "Updated raw data to actually include all samples. Specifically, Spleen_D30 samples from 527_EU_Suthar_60RNAseq..."
)


# ---- set up data frame ----
# Identify count columns

expr <- load_checkpoint(
  "all_samples_updated",
  dir = "projects/mavs_mito/data/r_objects"
)


rename_map <- c(
  "6_X_Spleen_D30" = "6_Spleen_D30",
  "18_Y_Spleen_D30" = "18_Spleen_D30",
  "30_Z_Spleen_D30" = "30_Spleen_D30",
  "42_A_Spleen_D30" = "42_Spleen_D30",
  "54_B_Spleen_D30" = "54_Spleen_D30"
  # etc.
)

# Rename only columns that exist
old <- colnames(expr)
hits <- intersect(old, names(rename_map))
colnames(expr)[match(hits, old)] <- unname(rename_map[hits])

# make counts object minus first column of expr
counts <- expr[, -1]
counts[is.na(counts)] <- 0 # replace NA values with 0
unique(is.na(counts$`1_WT_d0`)) # check for NA values in the renamed column

# ---- QC For samples ----
qc_list <- lapply(colnames(counts), function(s) {
  x <- counts[, s, drop = TRUE]

  list(
    qc = data.frame(
      sample = s,
      total_reads = sum(x),
      detected_genes = sum(x > 0),
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

list(
  qc_table = qc_summary,
  gene_summaries = gene_summaries
)

# Drop low-count samples based on qc_summary/qc_list above
counts <- counts[, colnames(counts) != "6_Spleen_D30"]
counts <- counts[, colnames(counts) != "42_Spleen_D30"]
counts <- counts[, colnames(counts) != "18_Spleen_D30"]

#drop old mock
counts <- counts[, colnames(counts) != "1_Spleen_D0"]
counts <- counts[, colnames(counts) != "13_Spleen_D0"]
counts <- counts[, colnames(counts) != "25_Spleen_D0"]
counts <- counts[, colnames(counts) != "37_Spleen_D0"]
counts <- counts[, colnames(counts) != "49_Spleen_D0"]

counts <- counts[, colnames(counts) != "10_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "22_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "34_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "46_SpleenMAVSko_D0"]
counts <- counts[, colnames(counts) != "58_SpleenMAVSko_D0"]

# Save/load checkpoint
save_checkpoint(counts, "counts_filtered", dir = "projects/mavs_mito/data/r_objects",
                lines = c(225:293),
                notes = "Filtered out Spleen samples 6, 18, and 42 at day 30 due to low counts and removed old mock"
)

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
                lines = c(308:343),
                notes = "Meta data for spleen D0 vs D7 samples"
)

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
  lines = c(309:370)
)

save_checkpoint(res_voom_spleen_d0_d7, "res_voom_spleen_d0_d7",
  notes = "limma-voom DE results for spleen D0 vs D7; named list of topTables (within-genotype day effects, direct genotype comparison at D7, interaction). Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(309:370)
)

save_checkpoint(res_voom_spleen_d0_d7_annotated, "res_voom_spleen_d0_d7_annotated",
  notes = "Gene Symbol annotated limma-voom DE results for spleen D0 vs D7. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(309:370)
)

# ---- GSEA Spleen ----
library(fgsea)

res_voom_spleen_d0_d7_annotated <- load_checkpoint("res_voom_spleen_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects")

gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

res_wt_vs_mavsko_d7 <- res_voom_spleen_d0_d7_annotated[["WT_vs_MAVSKO_D7"]]

fgsea_wt_vs_mavsko_d7 <- run_gsea(res_wt_vs_mavsko_d7, pathways)

# sanity checks
fgsea_wt_vs_mavsko_d7 %>% select(pathway, NES, pval, padj) %>% head(15)

fgsea_wt_vs_mavsko_d7 %>% filter(grepl("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|GLYCOLYSIS|HYPOXIA", pathway)) %>%
  select(pathway, NES, pval, padj)

sum(fgsea_wt_vs_mavsko_d7$padj < 0.05)

save_checkpoint(fgsea_wt_vs_mavsko_d7, "fgsea_wt_vs_mavsko_d7",
  notes = "GSEA results for WT vs MAVSKO at D7 in spleen. Pathways from MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(394:412)
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
                lines = c(463:496),
                notes = "Meta data for brain D0 vs D7 samples"
)

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
  lines = c(463:549)
)

save_checkpoint(res_voom_brain_d0_d7, "res_voom_brain_d0_d7",
  notes = "limma-voom DE results for brain D0 vs D7; named list of topTables (within-genotype day effects, direct genotype comparison at D7, interaction). Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(463:549)
)

save_checkpoint(res_voom_brain_d0_d7_annotated, "res_voom_brain_d0_d7_annotated",
  notes = "Gene Symbol annotated limma-voom DE results for brain D0 vs D7. Default filtering 70% of samples in any group. min count = 10",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(463:549)
)

# ---- GSEA ----
res_voom_brain_d0_d7_annotated <- load_checkpoint("res_voom_brain_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects")

gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

res_wt_vs_mavsko_d7_brain <- res_voom_brain_d0_d7_annotated[["WT_vs_MAVSKO_D7"]]

fgsea_wt_vs_mavsko_d7_brain <- run_gsea(res_wt_vs_mavsko_d7_brain, pathways)

# sanity checks
fgsea_wt_vs_mavsko_d7_brain %>% select(pathway, NES, pval, padj) %>% head(15)

fgsea_wt_vs_mavsko_d7_brain %>% filter(grepl("OXIDATIVE_PHOSPHORYLATION|FATTY_ACID|GLYCOLYSIS|HYPOXIA", pathway)) %>%
  select(pathway, NES, pval, padj)

sum(fgsea_wt_vs_mavsko_d7_brain$padj < 0.05)

save_checkpoint(fgsea_wt_vs_mavsko_d7_brain, "fgsea_wt_vs_mavsko_d7_brain",
  notes = "GSEA results for WT vs MAVSKO at D7 in brain. Pathways from MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(551:572)
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
                lines = c(614:636),
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
  lines = c(614:672)
)

save_checkpoint(res_voom_d7_tissue, "res_voom_d7_tissue",
  notes = "limma-voom DE: Brain - Spleen at D7 within WT and within MAVSKO (Spleen = reference). Named list of topTables.",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(614:672)
)

save_checkpoint(res_voom_d7_tissue_annotated, "res_voom_d7_tissue_annotated",
  notes = "Gene Symbol annotated Brain - Spleen D7 DE (WT and MAVSKO).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(614:672)
)

# ---- GSEA: Brain vs Spleen D7 (Hallmark) ----
res_voom_d7_tissue_annotated <- load_checkpoint("res_voom_d7_tissue_annotated",
  dir = "projects/mavs_mito/data/r_objects")

gmt_path <- "genesets/mh.all.v2026.1.Mm.symbols.gmt"
pathways <- gmtPathways(gmt_path)

fgsea_wt_brain_vs_spleen_d7     <- run_gsea(res_voom_d7_tissue_annotated[["WT_Brain_vs_Spleen"]], pathways)
fgsea_mavsko_brain_vs_spleen_d7 <- run_gsea(res_voom_d7_tissue_annotated[["MAVSKO_Brain_vs_Spleen"]], pathways)

# sanity checks (positive NES = enriched in Brain relative to Spleen)
fgsea_wt_brain_vs_spleen_d7 %>% select(pathway, NES, pval, padj) %>% head(15)
fgsea_mavsko_brain_vs_spleen_d7 %>% select(pathway, NES, pval, padj) %>% head(15)

sum(fgsea_wt_brain_vs_spleen_d7$padj < 0.05)
sum(fgsea_mavsko_brain_vs_spleen_d7$padj < 0.05)

save_checkpoint(fgsea_wt_brain_vs_spleen_d7, "fgsea_wt_brain_vs_spleen_d7",
  notes = "GSEA Brain vs Spleen at D7, WT (Spleen = ref). MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(688:703)
)

save_checkpoint(fgsea_mavsko_brain_vs_spleen_d7, "fgsea_mavsko_brain_vs_spleen_d7",
  notes = "GSEA Brain vs Spleen at D7, MAVSKO (Spleen = ref). MSigDB Hallmark gene sets (v2026.1).",
  dir = "projects/mavs_mito/data/r_objects",
  lines = c(688:703)
)

# load checkpoint
