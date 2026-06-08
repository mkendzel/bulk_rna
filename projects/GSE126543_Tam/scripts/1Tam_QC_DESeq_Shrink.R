# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)
library(ggvenn)
library(grid)
library(GEOquery)
library(readxl)
# ---- Data ----


# Load metadata
gse <- getGEO(filename = "projects/GSE126543_Tam/data/GSE124439_series_matrix.txt")
metadata <- pData(gse)
colnames(metadata)
clinical <- read_excel("projects/GSE126543_Tam/data/1-s2.0-S221112471931263X-mmc3.xlsx",
                       sheet = "Table S1A")

head(clinical)
colnames(clinical)
# Load counts
counts <- read.table("projects/GSE126543_Tam/data/GSE124439_raw_counts_GRCh38.p13_NCBI.tsv",
                     header = TRUE,
                     sep = "\t",
                     row.names = 1,
                     check.names = FALSE)

dim(counts)
head(counts[, 1:5])

# Check subregion by group together
table(metadata$`cns subregion:ch1`, metadata$`sample group:ch1`)

# ---- Subset to only motor cortex samples ----
# Lateral Motor Cortex
meta_lateral <- metadata[metadata$`sample group:ch1` %in% c("ALS Spectrum MND", "Non-Neurological Control") &
                           metadata$`cns subregion:ch1` == "Motor Cortex (Lateral)", ]

counts_lateral <- counts[, colnames(counts) %in% rownames(meta_lateral)]

# Medial Motor Cortex
meta_medial <- metadata[metadata$`sample group:ch1` %in% c("ALS Spectrum MND", "Non-Neurological Control") &
                          metadata$`cns subregion:ch1` == "Motor Cortex (Medial)", ]

counts_medial <- counts[, colnames(counts) %in% rownames(meta_medial)]

# Check sample numbers for each
table(meta_lateral$`sample group:ch1`)
table(meta_medial$`sample group:ch1`)

# Check dimensions
dim(counts_lateral)
dim(counts_medial)



clinical$`RNA-seq ID` <- trimws(clinical$`RNA-seq ID`)
meta_lateral$title <- trimws(meta_lateral$title)
meta_medial$title <- trimws(meta_medial$title)

# Merge on title = RNA-seq ID
meta_lateral <- merge(meta_lateral,
                      clinical[, c("RNA-seq ID", "Gender", "ALS NMF subtype**")],
                      by.x = "title",
                      by.y = "RNA-seq ID",
                      all.x = TRUE)

meta_medial <- merge(meta_medial,
                     clinical[, c("RNA-seq ID", "Gender", "ALS NMF subtype**")],
                     by.x = "title",
                     by.y = "RNA-seq ID",
                     all.x = TRUE)

table(meta_lateral$Gender)
table(meta_lateral$`ALS NMF subtype**`)
table(meta_medial$Gender)
table(meta_medial$`ALS NMF subtype**`)


# ---- QC For samples ----

#Save/load checkpoint
# saveRDS(counts, "projects/GSE126543_Tam/data/r_objects/counts_subset(motorcortex).rds")
# counts <- readRDS("projects/GSE126543_Tam/data/r_objects/counts_subset(motorcortex).rds")

# ---- Lateral Motor Cortex dds----
meta_lateral <- meta_lateral[match(colnames(counts_lateral), meta_lateral$geo_accession), ]
all(meta_lateral$geo_accession == colnames(counts_lateral))

meta_lateral$condition <- ifelse(meta_lateral$`sample group:ch1` == "ALS Spectrum MND",
                                 "ALS", "Control")
meta_lateral$condition <- factor(meta_lateral$condition, levels = c("Control", "ALS"))
meta_lateral$Gender <- factor(meta_lateral$Gender)

table(meta_lateral$condition)
table(meta_lateral$Gender)

dds_lateral <- DESeqDataSetFromMatrix(countData = counts_lateral,
                                      colData   = meta_lateral,
                                      design    = ~ Gender + condition)

# ---- Medial Motor Cortex dds----
meta_medial <- meta_medial[match(colnames(counts_medial), meta_medial$geo_accession), ]
all(meta_medial$geo_accession == colnames(counts_medial))

meta_medial$condition <- ifelse(meta_medial$`sample group:ch1` == "ALS Spectrum MND",
                                "ALS", "Control")
meta_medial$condition <- factor(meta_medial$condition, levels = c("Control", "ALS"))
meta_medial$Gender <- factor(meta_medial$Gender)

table(meta_medial$condition)
table(meta_medial$Gender)

dds_medial <- DESeqDataSetFromMatrix(countData = counts_medial,
                                     colData   = meta_medial,
                                     design    = ~ Gender + condition)

# ---- Lateral Motor Cortex Deseq ----
keep_lateral <- rowSums(DESeq2::counts(dds_lateral) >= 10) >= 5
dds_lateral <- dds_lateral[keep_lateral, ]
vsd_lateral <- DESeq2::vst(dds_lateral)
DESeq2::plotPCA(vsd_lateral, intgroup = c("condition", "Gender")) + 
  ggplot2::ggtitle("Motor Cortex (Lateral)")

dds_lateral <- DESeq2::DESeq(dds_lateral)

res_lateral <- DESeq2::results(dds_lateral)
# Order by adjusted p-value
res_lateral <- res_lateral[order(res_lateral$padj), ]
head(res_lateral, 20)
summary(res_lateral)

# ---- Medial Motor Cortex DESeq ----
keep_medial <- rowSums(DESeq2::counts(dds_medial) >= 10) >= 5
dds_medial <- dds_medial[keep_medial, ]
vsd_medial <- DESeq2::vst(dds_medial)
DESeq2::plotPCA(vsd_medial, intgroup = c("condition", "Gender")) + 
  ggplot2::ggtitle("Motor Cortex (Medial)")

dds_medial <- DESeq2::DESeq(dds_medial)
results(dds_medial)

# Fit the DESeq2 model; disable outlier replacement for consistent behavior across conditions
dds <- DESeq2::DESeq(dds, minReplicatesForReplace = Inf)

# Inspect the coefficient names to confirm which contrasts exist
coef_names <- DESeq2::resultsNames(dds)

# Keep only coefficients that represent contrasts against the Mock_0 reference
coef_names_vs_mock <- coef_names[grep("vs_Mock", coef_names)]

# Shrink log2 fold-changes with apeglm for each vs-Mock coefficient
res_shrunk_list <- setNames(
  lapply(coef_names_vs_mock, function(coef_name) {
    DESeq2::lfcShrink(dds, coef = coef_name, type = "apeglm")
  }),
  coef_names_vs_mock
)

# Quick check of unshrunk MA plot for the default results
res <- DESeq2::results(dds)
DESeq2::plotMA(res, ylim = c(-5, 5))


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
# ---- Collect all Ensembl IDs across contrasts ----
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














