# ---- Libraries ----
library(dplyr)
library(readxl)
library(limma)
library(edgeR)
library(stringr)
library(ggvenn)
library(grid)


# ----- Import data ----
## Shilu's data
#Import expression matrix in tsv format
expr_R2SDHF <- read.delim(
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
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


# Build metadata
meta_R2SDHF <- data.frame(
  sample = colnames(counts_R2SDHF),
  treatment = gsub("_\\d+_\\d+$", "", colnames(counts_R2SDHF)),
  time = as.numeric(gsub("^[A-Za-z]+_(\\d+)_\\d+$", "\\1", colnames(counts_R2SDHF))),
  rep = gsub("^.*_(\\d+)$", "\\1", colnames(counts_R2SDHF))
)
rownames(meta_R2SDHF) <- meta_R2SDHF$sample

# Process
dge_R2SDHF <- DGEList(counts = counts_R2SDHF)
keep_R2SDHF <- filterByExpr(dge_R2SDHF)
dge_R2SDHF <- dge_R2SDHF[keep_R2SDHF, , keep.lib.sizes = FALSE]
dge_R2SDHF <- calcNormFactors(dge_R2SDHF)

# Simple design for voom
design_R2SDHF <- model.matrix(~ treatment * factor(time), data = meta_R2SDHF)
v_R2SDHF <- voom(dge_R2SDHF, design_R2SDHF)

# ---- Combined PCA ----
# Find shared genes
shared_genes <- intersect(rownames(v$E), rownames(v_R2SDHF$E))

# Combine expression
combined_expr <- cbind(v$E[shared_genes, ], v_R2SDHF$E[shared_genes, ])

# Batch correct
batch <- factor(c(rep("organoid", ncol(v$E)), rep("R2SDHF", ncol(v_R2SDHF$E))))
corrected <- removeBatchEffect(combined_expr, batch = batch)

# PCA
pca_combined <- prcomp(t(corrected), scale. = TRUE)

# Build plotting df
pca_combined_df <- data.frame(
  PC1 = pca_combined$x[, 1],
  PC2 = pca_combined$x[, 2],
  dataset = c(rep("organoid", ncol(v$E)), rep("R2SDHF", ncol(v_R2SDHF$E))),
  label = c(as.character(meta_filtered$Stage), paste0(meta_R2SDHF$treatment, "_", meta_R2SDHF$time, "h"))
)

ggplot(pca_combined_df, aes(x = PC1, y = PC2, color = label, shape = dataset)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", round(summary(pca_combined)$importance[2,1]*100, 1), "%)"),
    y = paste0("PC2 (", round(summary(pca_combined)$importance[2,2]*100, 1), "%)")
  )
##### Jimina's data
#### Expression matrix
# Define base path
base_path <- "data/Anderson-Suthar_collab/Als_org_Jimena"

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
  select(Admera_Health_ID, Patient, AgeOrg, Line, Replicate, Sex)

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

meta_filtered$Stage <- case_when(
  meta_filtered$AgeOrg %in% c(30,31,34,35) ~ "d30",
  meta_filtered$AgeOrg %in% c(50,52) ~ "d50",
  meta_filtered$AgeOrg %in% c(75,76,77) ~ "d75",
  meta_filtered$AgeOrg == 120 ~ "d120"
)

meta_filtered$Stage <- factor(meta_filtered$Stage, levels = c("d30", "d50", "d75", "d120"))

#Count matrix
counts <- expr %>% select(all_of(meta_filtered$Admera_Health_ID))

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
saveRDS(counts, "data/Anderson-Suthar_collab/Als_org_Jimena/counts(hSpS_Ctrx).rds")
counts <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/counts(hSpS_Ctrx).rds")

saveRDS(meta_filtered, "data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx).rds")

table(meta_filtered$Stage, meta_filtered$Line)
# ---- limma and edgeR set up ----

dge <- DGEList(counts = counts)

# Filter low expressed genes - keep genes expressed in at least a reasonable number of samples
keep <- filterByExpr(dge, group = meta_filtered$Line)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Set up design matrix
design <- model.matrix(~ Sex + Stage * Line, data = meta_filtered)
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

# ---- View Results ----
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

#Save/load DESeq object
# saveRDS(dds, "data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")
dds <- readRDS("data/R2SDHF/r_objects/plasmo_dds_noMAD.rds")

# Normalied reads
dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)


#Save/load Shrunk Data
saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo_resShrink_noMAD_qcmin10.rds")
res_shrunk_list <- readRDS("data/R2SDHF/Plasmo_resShrink_noMAD_qcmin10.rds")


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

saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")
res_shrunk_list <- readRDS("data/R2SDHF/Plasmo/Plasmo_resShrink_noMAD_qcmin10_annotated.rds")














