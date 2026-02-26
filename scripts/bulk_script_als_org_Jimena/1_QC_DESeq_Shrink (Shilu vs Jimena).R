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
expr <- read.delim(
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
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

meta_filtered$Line <- relevel(meta_filtered$Line, ref = "control")
meta_filtered$AgeOrg <- as.factor(meta_filtered$AgeOrg)
meta_filtered$Patient <- as.factor(meta_filtered$Patient)
meta_filtered$Sex <- as.factor(meta_filtered$Sex)

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
# saveRDS(counts, "data/Anderson-Suthar_collab/Als_org_Jimena/counts(hSpS_Ctrx).rds")
# counts <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/counts(hSpS_Ctrx).rds")

saveRDS(meta_filtered, "data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx).rds")
meta_filtered <- readRDS("data/Anderson-Suthar_collab/Als_org_Jimena/meta_filtered(hSpS_Ctrx).rds")

# ---- limma and edgeR set up ----

dge <- DGEList(counts = counts)

# Filter low expressed genes - keep genes expressed in at least a reasonable number of samples
keep <- filterByExpr(dge, group = meta_filtered$Line)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# Normalize
dge <- calcNormFactors(dge)

# Set up design matrix
design <- model.matrix(~ Sex + AgeOrg + Line, data = meta_filtered)

# First voom transformation
v <- voom(dge, design)

# Account for repeated measures (Patient)
corfit <- duplicateCorrelation(v, design, block = meta_filtered$Patient)

# Second voom with the correlation estimate
v <- voom(dge, design, block = meta_filtered$Patient, correlation = corfit$consensus)

# Re-estimate correlation with updated voom weights
corfit <- duplicateCorrelation(v, design, block = meta_filtered$Patient)

# Fit the model
fit <- lmFit(v, design, block = meta_filtered$Patient, correlation = corfit$consensus)
fit <- eBayes(fit)

# Get results for ALS vs Control
results <- topTable(fit, coef = "LineALS", number = Inf)










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














