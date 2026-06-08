library(tidyverse)
library(edgeR)
library(limma)
library(org.Hs.eg.db)
library(AnnotationDbi)
# Import expression matrix
expr_R2SDHF <- read.delim(
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Identify count columns
count_cols_R2SDHF <- grep("_count$", colnames(expr_R2SDHF), value = TRUE)
counts_R2SDHF <- expr_R2SDHF[, c("gene_id", count_cols_R2SDHF)]
rownames(counts_R2SDHF) <- counts_R2SDHF$gene_id
counts_R2SDHF$gene_id <- NULL
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

# Remove MAD only, keep ROV
counts_R2SDHF <- counts_R2SDHF[, !grepl("^MAD_", colnames(counts_R2SDHF))]

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

# ---- LIMMA voom analysis ----
design <- model.matrix(~ 0 + Group, data = meta_R2SDHF)
colnames(design) <- levels(meta_R2SDHF$Group)

dge <- DGEList(counts = counts_R2SDHF)
keep <- filterByExpr(dge, design)
dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- calcNormFactors(dge)

v <- voom(dge, design)
fit <- lmFit(v, design)

# Define contrasts including ROV
contr <- makeContrasts(
  TX12_vs_Mock  = TX_12 - Mock_0,
  TX24_vs_Mock  = TX_24 - Mock_0,
  TX48_vs_Mock  = TX_48 - Mock_0,
  TX48_vs_TX12  = TX_48 - TX_12,
  ROV12_vs_Mock = ROV_12 - Mock_0,
  ROV24_vs_Mock = ROV_24 - Mock_0,
  ROV48_vs_Mock = ROV_48 - Mock_0,
  levels = design
)

fit2 <- contrasts.fit(fit, contr)
fit2 <- eBayes(fit2)


# ROV results
res_ROV12 <- topTable(fit2, coef = "ROV12_vs_Mock", number = Inf)
res_ROV24 <- topTable(fit2, coef = "ROV24_vs_Mock", number = Inf)
res_ROV48 <- topTable(fit2, coef = "ROV48_vs_Mock", number = Inf)

saveRDS(fit2, "data/Anderson-Suthar_collab/Als_org_Jimena/fit2(R2SDHF).rds")
saveRDS(v, "data/Anderson-Suthar_collab/Als_org_Jimena/voom(R2SDHF).rds")
saveRDS(res_ROV12, "data/Anderson-Suthar_collab/r_objects/res_ROV12(R2SDHF).rds")
saveRDS(res_ROV24, "data/Anderson-Suthar_collab/r_objects/res_ROV24(R2SDHF).rds")
saveRDS(res_ROV48, "data/Anderson-Suthar_collab/r_objects/res_ROV48(R2SDHF).rds")