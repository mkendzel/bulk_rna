# ---- Libraries ----
library(dplyr)
library(readxl)
library(limma)
library(edgeR)
library(stringr)
library(ggvenn)
library(grid)
library(org.Hs.eg.db)
library(uwot)

#import all helper functions in /R
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# ----- Import and process R2SDHF data ----
## Shilu's data
#Import expression matrix in tsv format
expr_R2SDHF <- read.delim(
  "projects/R2SDHF/data/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify count columns
count_cols_R2SDHF <- grep("_count$", colnames(expr_R2SDHF), value = TRUE)
# Subset
counts_R2SDHF <- expr_R2SDHF[, c("gene_id", count_cols_R2SDHF)]
# Set rownames
rownames(counts_R2SDHF) <- counts_R2SDHF$gene_id
counts_R2SDHF$gene_id <- NULL
# Remove "_count" suffix
colnames(counts_R2SDHF) <- sub("_count$", "", colnames(counts_R2SDHF))
# Rename columns
rename_map_R2SDHF <- c(
  "R2SDHF_1"  = "Mock_0_1",
  "R2SDHF_2"  = "Mock_0_2",
  "R2SDHF_3"  = "Mock_0_3",
  "R2SDHF_4"  = "TX_12_1",
  "R2SDHF_5"  = "TX_12_2",
  "R2SDHF_6"  = "TX_12_3",
  "R2SDHF_7"  = "ROV_12_1",
  "R2SDHF_8"  = "ROV_12_2",
  "R2SDHF_9"  = "ROV_12_3",
  "R2SDHF_10" = "MAD_12_1",
  "R2SDHF_11" = "MAD_12_2",
  "R2SDHF_12" = "MAD_12_3",
  "R2SDHF_13" = "TX_24_1",
  "R2SDHF_14" = "TX_24_2",
  "R2SDHF_15" = "TX_24_3",
  "R2SDHF_16" = "ROV_24_1",
  "R2SDHF_17" = "ROV_24_2",
  "R2SDHF_18" = "ROV_24_3",
  "R2SDHF_19" = "MAD_24_1",
  "R2SDHF_20" = "MAD_24_2",
  "R2SDHF_21" = "MAD_24_3",
  "R2SDHF_22" = "TX_48_1",
  "R2SDHF_23" = "TX_48_2",
  "R2SDHF_24" = "TX_48_3",
  "R2SDHF_25" = "ROV_48_1",
  "R2SDHF_26" = "ROV_48_2",
  "R2SDHF_27" = "ROV_48_3",
  "R2SDHF_28" = "MAD_48_1",
  "R2SDHF_29" = "MAD_48_2",
  "R2SDHF_30" = "MAD_48_3"
)
colnames(counts_R2SDHF) <- rename_map_R2SDHF[colnames(counts_R2SDHF)]

counts_R2SDHF <- counts_R2SDHF[, !grepl("^MAD_", colnames(counts_R2SDHF))]
# counts_R2SDHF <- counts_R2SDHF[, !grepl("^ROV_", colnames(counts_R2SDHF))]

# Convert Ensembl IDs to gene symbols
gene_map <- mapIds(org.Hs.eg.db,
                   keys = rownames(counts_R2SDHF),
                   keytype = "ENSEMBL",
                   column = "SYMBOL",
                   multiVals = "first")

cat("Unmapped:", sum(is.na(gene_map)), "\n")
cat("Duplicate symbols:", sum(duplicated(gene_map[!is.na(gene_map)])), "\n")

# Keep highest-expressed Ensembl ID per symbol
mapped <- data.frame(
  ensembl = names(gene_map),
  symbol = gene_map,
  mean_expr = rowMeans(counts_R2SDHF)
) %>%
  filter(!is.na(symbol)) %>%
  arrange(symbol, desc(mean_expr)) %>%
  filter(!duplicated(symbol))

counts_R2SDHF <- counts_R2SDHF[mapped$ensembl, ]
rownames(counts_R2SDHF) <- mapped$symbol

# Build metadata
meta_R2SDHF <- data.frame(
  sample = colnames(counts_R2SDHF),
  treatment = gsub("_\\d+_\\d+$", "", colnames(counts_R2SDHF)),
  time = gsub("^[A-Za-z]+_(\\d+)_\\d+$", "\\1", colnames(counts_R2SDHF)),
  rep = gsub("^.*_(\\d+)$", "\\1", colnames(counts_R2SDHF))
)
rownames(meta_R2SDHF) <- meta_R2SDHF$sample

meta_R2SDHF$treatment <- factor(meta_R2SDHF$treatment, levels = c("Mock", "TX", "ROV"))
meta_R2SDHF$time <- factor(meta_R2SDHF$time)
meta_R2SDHF$Group <- factor(paste0(meta_R2SDHF$treatment, "_", meta_R2SDHF$time))

# Combined condition with Mock_0 as reference
meta_R2SDHF$condition <- ifelse(meta_R2SDHF$treatment == "Mock", "Mock_0",
                                paste0(meta_R2SDHF$treatment, "_", meta_R2SDHF$time))
meta_R2SDHF$condition <- relevel(factor(meta_R2SDHF$condition), ref = "Mock_0")

# ---- R2SDHF DE analysis ----
dge_R2SDHF <- DGEList(counts = counts_R2SDHF)
keep_R2SDHF <- filterByExpr(dge_R2SDHF, group = meta_R2SDHF$condition)
dge_R2SDHF <- dge_R2SDHF[keep_R2SDHF, , keep.lib.sizes = FALSE]
dge_R2SDHF <- calcNormFactors(dge_R2SDHF)

design_R2SDHF <- model.matrix(~ condition, data = meta_R2SDHF)
colnames(design_R2SDHF) <- make.names(colnames(design_R2SDHF))

v_R2SDHF <- voom(dge_R2SDHF, design_R2SDHF)
fit_R2SDHF <- lmFit(v_R2SDHF, design_R2SDHF)
fit_R2SDHF <- eBayes(fit_R2SDHF)

colnames(fit_R2SDHF$coefficients)

# map gene names
ensembl_ids <- rownames(v_R2SDHF$E)
gene_map <- mapIds(org.Hs.eg.db, keys = ensembl_ids, 
                   keytype = "ENSEMBL", column = "SYMBOL")

# Remove NAs and duplicates
mapped <- gene_map[!is.na(gene_map)]
mapped <- mapped[!duplicated(mapped)]

# Subset and rename
v_R2SDHF_sym <- v_R2SDHF$E[names(mapped), ]
rownames(v_R2SDHF_sym) <- mapped

### Save/load checkpoint
save_checkpoint(v_R2SDHF_sym, "voom_E_R2SDHF_sym",
                notes = "Voom expression matrix for R2SDHF with gene symbols as rownames after Ensembl mapping",
                dir = "projects/R2SDHF/data/r_objects", lines = 128:138)
# v_R2SDHF_sym <- load_checkpoint("voom_E_R2SDHF_sym", dir = "projects/R2SDHF/data/r_objects")

# ---- R2SDHF LIMMA voom analysis ----
design <- model.matrix(~ 0 + Group, data = meta_R2SDHF)
colnames(design) <- levels(meta_R2SDHF$Group)

#DGE list
dge <- DGEList(counts = counts_R2SDHF)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)
# Voom normalization (no blocking needed)
v <- voom(dge, design)

# Fit
fit <- lmFit(v, design)

# Define contrasts of interest
contr <- makeContrasts(
  TX12_vs_Mock = TX_12 - Mock_0,
  TX24_vs_Mock = TX_24 - Mock_0,
  TX48_vs_Mock = TX_48 - Mock_0,
  ROV12_vs_Mock = ROV_12 - Mock_0,
  ROV24_vs_Mock = ROV_24 - Mock_0,
  ROV48_vs_Mock = ROV_48 - Mock_0,
  ROV48_vs_ROV12 = ROV_48 - ROV_12,
  TX48_vs_TX12 = TX_48 - TX_12,
  levels = design
)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)

# Results for one contrast at a time
res_TX12 <- topTable(fit2, coef = "TX12_vs_Mock", number = Inf)
res_TX24 <- topTable(fit2, coef = "TX24_vs_Mock", number = Inf)
res_TX48 <- topTable(fit2, coef = "TX48_vs_Mock", number = Inf)
res_ROV12 <- topTable(fit2, coef = "ROV12_vs_Mock", number = Inf)
res_ROV24 <- topTable(fit2, coef = "ROV24_vs_Mock", number = Inf)
res_ROV48 <- topTable(fit2, coef = "ROV48_vs_Mock", number = Inf)
res_ROV48_vs_ROV12 <- topTable(fit2, coef = "ROV48_vs_ROV12", number = Inf)
res_TX48_vs_TX12 <- topTable(fit2, coef = "TX48_vs_TX12", number = Inf)

save_checkpoint(fit2, "fit2_R2SDHF",
                notes = "Limma fit object for R2SDHF data with contrasts for each treatment vs Mock and some pairwise contrasts",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:173))
save_checkpoint(v, "voom_R2SDHF",
                notes = "Voom object for R2SDHF data after filtering and normalization, with gene symbols as rownames",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:173))
save_checkpoint(res_TX12, "res_TX12_R2SDHF",
                notes = "DE results for TX 12h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_TX24, "res_TX24_R2SDHF",
                notes = "DE results for TX 24h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_TX48, "res_TX48_R2SDHF",
                notes = "DE results for TX 48h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_ROV12, "res_ROV12_R2SDHF",
                notes = "DE results for ROV 12h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_ROV24, "res_ROV24_R2SDHF",
                notes = "DE results for ROV 24h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_ROV48, "res_ROV48_R2SDHF",
                notes = "DE results for ROV 48h vs Mock in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_ROV48_vs_ROV12, "res_ROV48_vs_ROV12_R2SDHF",
                notes = "DE results for ROV 48h vs ROV 12h in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))
save_checkpoint(res_TX48_vs_TX12, "res_TX48_vs_TX12_R2SDHF",
                notes = "DE results for TX 48h vs TX 12h in R2SDHF data from limma fit",
                dir = "projects/R2SDHF/data/r_objects", lines = c(145:183))

# ---- Reload checkpoints (optional) ----
# fit2          <- load_checkpoint("fit2_R2SDHF",            dir = "projects/R2SDHF/data/r_objects")
# v             <- load_checkpoint("voom_R2SDHF",            dir = "projects/R2SDHF/data/r_objects")
# res_TX12      <- load_checkpoint("res_TX12_R2SDHF",        dir = "projects/R2SDHF/data/r_objects")
# res_TX24      <- load_checkpoint("res_TX24_R2SDHF",        dir = "projects/R2SDHF/data/r_objects")
# res_TX48      <- load_checkpoint("res_TX48_R2SDHF",        dir = "projects/R2SDHF/data/r_objects")
# res_ROV12     <- load_checkpoint("res_ROV12_R2SDHF",       dir = "projects/R2SDHF/data/r_objects")
# res_ROV24     <- load_checkpoint("res_ROV24_R2SDHF",       dir = "projects/R2SDHF/data/r_objects")
# res_ROV48     <- load_checkpoint("res_ROV48_R2SDHF",       dir = "projects/R2SDHF/data/r_objects")
# res_ROV48_vs_ROV12 <- load_checkpoint("res_ROV48_vs_ROV12_R2SDHF", dir = "projects/R2SDHF/data/r_objects")
# res_TX48_vs_TX12   <- load_checkpoint("res_TX48_vs_TX12_R2SDHF",   dir = "projects/R2SDHF/data/r_objects")


# ---- FgSEA for R2SDHF ----
library(fgsea)
library(msigdbr)

# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
pathways <- gmtPathways(gmt_path)


# Helper to build ranked list and run fgsea
run_gsea <- function(tt, pathways, nPermSimple = 10000) {
  # Resolve duplicates by keeping the entry with highest average expression
  n_before <- nrow(tt)
  tt <- tt[order(-tt$AveExpr), ]
  tt <- tt[!duplicated(rownames(tt)), ]
  n_dupes <- n_before - nrow(tt)
  message("Removed ", n_dupes, " duplicate gene names (kept highest AveExpr)")
  
  ranks <- setNames(tt$t, rownames(tt))
  ranks <- sort(ranks, decreasing = TRUE)
  
  res <- fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
  res <- res[order(res$pval), ]
  return(res)
}

fgsea_TX12 <- run_gsea(res_TX12, pathways)
fgsea_TX24 <- run_gsea(res_TX24, pathways)
fgsea_TX48 <- run_gsea(res_TX48, pathways)
fgsea_TX48_vs_TX12 <- run_gsea(res_TX48_vs_TX12, pathways)

save_checkpoint(fgsea_TX12, "fgsea_TX12_vs_Mock_R2SDHF",
                notes = "fgsea results for TX 12h vs Mock in R2SDHF data (Hallmark), ranked by limma t-statistic",
                dir = "projects/R2SDHF/data/r_objects", lines = 226:246)
save_checkpoint(fgsea_TX24, "fgsea_TX24_vs_Mock_R2SDHF",
                notes = "fgsea results for TX 24h vs Mock in R2SDHF data (Hallmark), ranked by limma t-statistic",
                dir = "projects/R2SDHF/data/r_objects", lines = 226:246)
save_checkpoint(fgsea_TX48, "fgsea_TX48_vs_Mock_R2SDHF",
                notes = "fgsea results for TX 48h vs Mock in R2SDHF data (Hallmark), ranked by limma t-statistic",
                dir = "projects/R2SDHF/data/r_objects", lines = 226:246)
save_checkpoint(fgsea_TX48_vs_TX12, "fgsea_TX48_vs_TX12_R2SDHF",
                notes = "fgsea results for TX 48h vs TX 12h in R2SDHF data (Hallmark), ranked by limma t-statistic",
                dir = "projects/R2SDHF/data/r_objects", lines = 226:246)
# ---- Jimina's data ----
#### Expression matrix
# Define base path
base_path <- "projects/Anderson-Suthar_collab/data/Als_org_Jimena"

# Import expression matrix
expr <- read_excel(
  file.path(base_path, "102025_lines1-6,c9_vs_ctrl,d30,50,75,120_allfeature_counts.xlsx")
)

# rowname set up
expr <- as.data.frame(expr)
rownames(expr) <- expr$Geneid
expr <- expr[, -1]

# colname formating
colnames(expr) <- sub("^\\./", "", colnames(expr))
colnames(expr) <- sub("_S.*$", "", colnames(expr))


#### Metadata
# Import metadata
metadata <- read_excel(
  file.path(base_path, "BulkSeqTargetList_tosend.xlsx")
)

metadata <- as.data.frame(metadata)

# ---- Set up and subset datasets ----

# Filter metadata for hSPS organoids with Ctrx coating
meta_filtered <- metadata %>%
  filter(Organoid == "hSpS", Coating == "Ctrx", Line %in% c("ALS", "control")) %>%
  dplyr::select(Admera_Health_ID, Patient, AgeOrg, Line, Replicate, Sex, Batch = d0)
meta_filtered %>%
  group_by(AgeOrg, Sex, Line) %>%
  summarise(
    n = n(),
    Patients = paste(unique(Patient), collapse = ", "),
    Lines = paste(unique(Line), collapse = ", ")
  )
meta_filtered$Line <- factor(meta_filtered$Line, levels = c("control", "ALS"))
meta_filtered$Line <- relevel(meta_filtered$Line, ref = "control")
meta_filtered$AgeOrg <- as.factor(meta_filtered$AgeOrg)
meta_filtered$Patient <- as.factor(meta_filtered$Patient)
meta_filtered$Sex <- as.factor(meta_filtered$Sex)
meta_filtered$Batch <- as.factor(meta_filtered$Batch)
meta_filtered$Stage <- case_when(
  meta_filtered$AgeOrg %in% c(30,31,34,35) ~ "d30",
  meta_filtered$AgeOrg %in% c(50,52) ~ "d50",
  meta_filtered$AgeOrg %in% c(75,76,77) ~ "d75",
  meta_filtered$AgeOrg == 120 ~ "d120"
)
meta_filtered$Stage <- factor(meta_filtered$Stage, levels = c("d30", "d50", "d75", "d120"))

#Count matrix
counts <- expr %>% dplyr::select(all_of(meta_filtered$Admera_Health_ID))

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

list(
  qc_table       = qc_summary,
  gene_summaries = gene_summaries
)

#Save/load checkpoint

save_checkpoint(counts, "filtered_counts_hSpS_Ctrx",
                notes = "Raw count matrix for hSPS organoids with Ctrx coating, subset to samples passing QC; rows=genes, cols=samples",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(220:275))

#optional load checkpoint
# counts <- load_checkpoint("filtered_counts_hSpS_Ctrx", dir = "projects/Anderson-Suthar_collab/data/r_objects")

save_checkpoint(meta_filtered, "meta_filtered_hSpS_Ctrx",
                notes = "Metadata for hSPS organoids with Ctrx coating, subset to samples passing QC; includes patient, line",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(241:272))
# optional load checkpoint
# meta_filtered <- load_checkpoint("meta_filtered_hSpS_Ctrx", dir = "projects/Anderson-Suthar_collab/data/r_objects")


# ---- limma and edgeR set up ----

dge <- DGEList(counts = counts)

# Filter low expressed genes - keep genes expressed in at least a reasonable number of samples
keep <- filterByExpr(dge, group = meta_filtered$Line)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Set u p design matrix
design <- model.matrix(~ Sex + Batch + Stage * Line, data = meta_filtered)
colnames(design) <- make.names(colnames(design))

# Voom + duplicateCorrelation (twice)
v <- voom(dge, design, block = meta_filtered$Patient)
corfit <- duplicateCorrelation(v, design, block = meta_filtered$Patient)
v <- voom(dge, design, block = meta_filtered$Patient, correlation = corfit$consensus)
corfit <- duplicateCorrelation(v, design, block = meta_filtered$Patient)

# Fit
fit <- lmFit(v, design, block = meta_filtered$Patient, correlation = corfit$consensus)
fit <- eBayes(fit)

# Check names
colnames(fit$coefficients)

# ---- Results ----

# Disease effect at each stage
con <- makeContrasts(
  d30_ALS = LineALS,
  d50_ALS = LineALS + Staged50.LineALS,
  d75_ALS = LineALS + Staged75.LineALS,
  d120_ALS = LineALS + Staged120.LineALS,
  levels = design
)
fit2 <- contrasts.fit(fit, con)
fit2 <- eBayes(fit2)

# Results per stage
res_d30 <- topTable(fit2, coef = "d30_ALS", number = Inf)
res_d50 <- topTable(fit2, coef = "d50_ALS", number = Inf)
res_d75 <- topTable(fit2, coef = "d75_ALS", number = Inf)
res_d120 <- topTable(fit2, coef = "d120_ALS", number = Inf)

# Save DE results
save_checkpoint(fit2, "fit2_hSpS_Ctrx_batch",
                notes = "Limma fit object with duplicateCorrelation for patient blocking, testing ALS vs control at each stage in hSPS organoids with Ctrx (coefs: d30_ALS, d50_ALS, d75_ALS, d120_ALS)",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:359))

save_checkpoint(v, "voom_hSpS_Ctrx_batch",
                notes = "Voom object for hSPS organoids with Ctrx coating, subset to samples passing QC",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:339))

save_checkpoint(res_d30, "res_d30_ALS_vs_Ctrl",
                notes = "DE results for ALS vs control at d30 in hSPS organoids with Ctrx coating, from limma fit with duplicateCorrelation for patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:365))
save_checkpoint(res_d50, "res_d50_ALS_vs_Ctrl",
                notes = "DE results for ALS vs control at d50 in hSPS organoids with Ctrx coating, from limma fit with duplicateCorrelation for patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:365))
save_checkpoint(res_d75, "res_d75_ALS_vs_Ctrl",
                notes = "DE results for ALS vs control at d75 in hSPS organoids with Ctrx coating, from limma fit with duplicateCorrelation for patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:365))
save_checkpoint(res_d120, "res_d120_ALS_vs_Ctrl",
                notes = "DE results for ALS vs control at d120 in hSPS organoids with Ctrx coating, from limma fit with duplicateCorrelation for patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(322:365))

# ---- Reload checkpoints (optional) ----
# fit2     <- load_checkpoint("fit2_hSpS_Ctrx_batch",  dir = "projects/Anderson-Suthar_collab/data/r_objects")
# v        <- load_checkpoint("voom_hSpS_Ctrx_batch",  dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d30  <- load_checkpoint("res_d30_ALS_vs_Ctrl",   dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d50  <- load_checkpoint("res_d50_ALS_vs_Ctrl",   dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d75  <- load_checkpoint("res_d75_ALS_vs_Ctrl",   dir = "projects/Anderson-Suthar_collab/data/r_objects")
res_d120 <- load_checkpoint("res_d120_ALS_vs_Ctrl",  dir = "projects/Anderson-Suthar_collab/data/r_objects")

# ---- GSEA ----
library(fgsea)
library(msigdbr)

# Get gene sets (example: Hallmark, human)
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
pathways <- gmtPathways(gmt_path)

# Add specific curated pathway(s) from MSigDB via msigdbr
# WP2447: ALS (Amyotrophic Lateral Sclerosis) - WikiPathways, C2 curated collection
wp2447_df <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:WIKIPATHWAYS") |>
  dplyr::filter(gs_exact_source == "WP2447")
pathways <- c(pathways, split(wp2447_df$gene_symbol, wp2447_df$gs_name))

# Helper to build ranked list and run fgsea
run_gsea <- function(tt, pathways, nPermSimple = 10000) {
  # Resolve duplicates by keeping the entry with highest average expression
  n_before <- nrow(tt)
  tt <- tt[order(-tt$AveExpr), ]
  tt <- tt[!duplicated(rownames(tt)), ]
  n_dupes <- n_before - nrow(tt)
  message("Removed ", n_dupes, " duplicate gene names (kept highest AveExpr)")
  
  ranks <- setNames(tt$t, rownames(tt))
  ranks <- sort(ranks, decreasing = TRUE)
  
  res <- fgsea(pathways = pathways, stats = ranks, nPermSimple = nPermSimple)
  res <- res[order(res$pval), ]
  return(res)
}

# Run for each timepoint
gsea_d30  <- run_gsea(res_d30, pathways)
gsea_d50  <- run_gsea(res_d50, pathways)
gsea_d75  <- run_gsea(res_d75, pathways)
gsea_d120 <- run_gsea(res_d120, pathways)

save_checkpoint(gsea_d30,  "fgsea_d30_ALS_vs_Ctrl",
                notes = "fgsea results for ALS vs control at d30 (Hallmark + WP2447), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(430:465))
save_checkpoint(gsea_d50,  "fgsea_d50_ALS_vs_Ctrl",
                notes = "fgsea results for ALS vs control at d50 (Hallmark + WP2447), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(430:465))
save_checkpoint(gsea_d75,  "fgsea_d75_ALS_vs_Ctrl",
                notes = "fgsea results for ALS vs control at d75 (Hallmark + WP2447), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(430:465))
save_checkpoint(gsea_d120, "fgsea_d120_ALS_vs_Ctrl",
                notes = "fgsea results for ALS vs control at d120 (Hallmark + WP2447), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = c(430:465))
# ---- GSEA networks ----
library(clusterProfiler)
library(enrichplot)
library(msigdbr)
library(org.Hs.eg.db)

# Get Hallmark gene sets formatted for clusterProfiler
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
pathways <- gmtPathways(gmt_path)

# Convert gmt pathways to t2g format
msig_t2g <- data.frame(
  gs_name = rep(names(pathways), lengths(pathways)),
  gene_symbol = unlist(pathways)
)

# Helper to build ranked list from topTable results
make_ranks <- function(tt) {
  ranks <- setNames(tt$t, rownames(tt))
  ranks <- sort(ranks, decreasing = TRUE)
  ranks <- ranks[!duplicated(names(ranks))]
  return(ranks)
}

# Run clusterProfiler GSEA for each timepoint
run_cp_gsea <- function(tt, t2g) {
  ranks <- make_ranks(tt)
  GSEA(geneList = ranks, TERM2GENE = t2g, pvalueCutoff = 1, eps = 0)
}

cp_d30  <- run_cp_gsea(res_d30, msig_t2g)
cp_d50  <- run_cp_gsea(res_d50, msig_t2g)
cp_d75  <- run_cp_gsea(res_d75, msig_t2g)
cp_d120 <- run_cp_gsea(res_d120, msig_t2g)

#Save
save_checkpoint(cp_d30,  "cp_gsea_d30_ALS_vs_Ctrl",
                notes = "clusterProfiler GSEA results for ALS vs control at d30 (Hallmark), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = 503:512)
save_checkpoint(cp_d50,  "cp_gsea_d50_ALS_vs_Ctrl",
                notes = "clusterProfiler GSEA results for ALS vs control at d50 (Hallmark), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = 503:512)
save_checkpoint(cp_d75,  "cp_gsea_d75_ALS_vs_Ctrl",
                notes = "clusterProfiler GSEA results for ALS vs control at d75 (Hallmark), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = 503:512)
save_checkpoint(cp_d120, "cp_gsea_d120_ALS_vs_Ctrl",
                notes = "clusterProfiler GSEA results for ALS vs control at d120 (Hallmark), ranked by limma t-statistic with patient blocking",
                dir = "projects/Anderson-Suthar_collab/data/r_objects", lines = 503:512)

# ---- graph Results ----
library(UpSetR)
sig_genes <- list(
  d30 = rownames(res_d30[res_d30$adj.P.Val < 0.05, ]),
  d50 = rownames(res_d50[res_d50$adj.P.Val < 0.05, ]),
  d75 = rownames(res_d75[res_d75$adj.P.Val < 0.05, ]),
  d120 = rownames(res_d120[res_d120$adj.P.Val < 0.05, ])
)
upset(fromList(sig_genes), order.by = "freq", 
      sets = c("d120", "d75", "d50", "d30"),
      keep.order = TRUE)

# PCA on the voom expression
pca <- prcomp(t(v$E), scale. = TRUE)

# Build a dataframe for plotting
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  Stage = meta_filtered$Stage,
  Line = meta_filtered$Line,
  Patient = meta_filtered$Patient
)

library(ggplot2)
ggplot(pca_df, aes(x = PC1, y = PC2, color = Stage, shape = Line)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca)$importance[2,2]*100, 1), "%)")
  )



# ---- Combined PCA ----
# Build biological group factor across both experiments
bio_group <- factor(c(
  ifelse(meta_filtered$Line == "ALS", "disease", "control"),
  ifelse(meta_R2SDHF$treatment == "Mock", "control", "disease")
))

design <- model.matrix(~bio_group)
corrected <- removeBatchEffect(combined_expr, batch = batch, design = design)
# PCA
pca_combined <- prcomp(t(corrected), scale. = TRUE)

# Build plotting df
pca_combined_df <- data.frame(
  PC1 = pca_combined$x[, 1],
  PC2 = pca_combined$x[, 2],
  dataset = c(rep("organoid", ncol(v$E)), rep("R2SDHF", ncol(v_R2SDHF_sym))),
  label = c(as.character(meta_filtered$Stage), paste0(meta_R2SDHF$treatment, "_", meta_R2SDHF$time, "h"))
)

library(RColorBrewer)

label_colors <- c(
  "d30"     = "#C6DBEF",
  "d50"     = "#6BAED6",
  "d75"     = "#2171B5",
  "d120"    = "#08306B",
  "Mock_0h" = "#2CA02C",
  "TX_12h"  = "#FCBBA1",
  "TX_24h"  = "#FB6A4A",
  "TX_48h"  = "#CB181D",
  "ROV_12h" = "#C994C7",
  "ROV_24h" = "#88419D",
  "ROV_48h" = "#4D004B"
)

ggplot(pca_combined_df, aes(x = PC1, y = PC2, color = label, shape = dataset)) +
  geom_point(size = 3) +
  scale_color_manual(values = label_colors) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_combined)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_combined)$importance[2,2]*100, 1), "%)")
  )


set.seed(42)
umap_combined <- umap(t(corrected), n_neighbors = 15, min_dist = 0.1)

umap_df <- data.frame(
  UMAP1 = umap_combined[, 1],
  UMAP2 = umap_combined[, 2],
  label = pca_combined_df$label,
  dataset = pca_combined_df$dataset
)

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = label, shape = dataset)) +
  geom_point(size = 3) +
  scale_color_manual(values = label_colors) +
  theme_minimal()


# =========================
# ---- RRHO2 ANalysis ----
# ========================
library(RRHO2)

# --- Experiment 1: Organoid (ALS vs control)
group_org <- factor(meta_filtered$Line, levels = c("control", "ALS"))
design_org <- model.matrix(~group_org)
fit_org <- lmFit(v, design_org)
fit_org <- eBayes(fit_org)
de_org <- topTable(fit_org, coef = 2, number = Inf)

# --- Experiment 2: R2SDHF (TX/ROV vs Mock) 
group_r2 <- factor(ifelse(meta_R2SDHF$treatment == "Mock", "Mock", "disease"), levels = c("Mock", "disease"))
design_r2 <- model.matrix(~group_r2)
fit_r2 <- lmFit(v_R2SDHF_sym, design_r2)
fit_r2 <- eBayes(fit_r2)
de_r2 <- topTable(fit_r2, coef = 2, number = Inf)

# RRHO2 expects a data frame with gene name and signed -log10(p-value)
shared <- intersect(rownames(de_org), rownames(de_r2))

list1 <- data.frame(
  gene = shared,
  value = -log10(de_org[shared, "P.Value"]) * sign(de_org[shared, "logFC"])
)

list2 <- data.frame(
  gene = shared,
  value = -log10(de_r2[shared, "P.Value"]) * sign(de_r2[shared, "logFC"])
)

rrho_res <- RRHO2_initialize(list1, list2, labels = c("Organoid ALS vs Ctrl", "R2SDHF Disease vs Mock"), log10.ind = TRUE)
RRHO2_heatmap(rrho_res)


# =========================
# ---- old ----
# ========================
 #Save/load DESeq object
# saveRDS(dds, "projects/R2SDHF/data/r_objects/plasmo_dds_noMAD.rds")
dds <- readRDS("projects/R2SDHF/data/r_objects/plasmo_dds_noMAD.rds")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#Save/load Shrunk Data
saveRDS(res_shrunk_list, "projects/R2SDHF/data/Plasmo_resShrink_noMAD_qcmin10.rds")
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo_resShrink_noMAD_qcmin10.rds")


# ---- ADD mapping of gene names ----
# Ensembl IDs across contrasts ----
all_ensembl <- res_shrunk_list |>
  purrr::map(~ rownames(.x)) |>
  unlist(use.names = FALSE) |>
  unique() |>
  as.character()

# Remove version suffix if present (e.g. ENSG00000123456.7)
all_ensembl_clean <- sub("\\..*$", "", all_ensembl)

# ---- Build mapping table ----
gene_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys    = all_ensembl_clean,
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
) |>
  as_tibble() |>
  distinct(ENSEMBL, .keep_all = TRUE) |>
  rename(
    ensembl_id = ENSEMBL,
    entrez_id  = ENTREZID,
    SYMBOL     = SYMBOL
  )

# ---- Append identifiers to every contrast ----
res_shrunk_list <- purrr::map(res_shrunk_list, function(df) {
  
  df |>
    as.data.frame() |>
    rownames_to_column("ensembl_id") |>
    mutate(ensembl_id = sub("\\..*$", "", ensembl_id)) |>
    left_join(gene_map, by = "ensembl_id") |>
    as_tibble()
  
})

saveRDS(res_shrunk_list, "projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("projects/R2SDHF/data/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")














