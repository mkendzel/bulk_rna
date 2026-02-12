# ---- Libraries ----
library(stringr)
library(dplyr)
library(DESeq2)

# ----- Import data ----
expr <- read.delim(
  "data/R2SDHF/Plasmo/R2SDHF-expression-matrix.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ---- set up data frame ----
# Identify count columns
count_cols <- grep("_count$", colnames(expr), value = TRUE)

# Subset
counts <- expr[, c("gene_id", count_cols)]

# Set rownames
rownames(counts) <- counts$gene_id
counts$gene_id <- NULL

# Remove "_count" suffix
colnames(counts) <- sub("_count$", "", colnames(counts))






# rename columns
# Map old -> new (explicit)
rename_map <- c(
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

# Rename only columns that exist
old <- colnames(counts)
hits <- intersect(old, names(rename_map))
colnames(counts)[match(hits, old)] <- unname(rename_map[hits])

# Checks
setdiff(names(rename_map), colnames(counts))   # should be character(0) if all expected cols present
any(duplicated(colnames(counts)))              # should be FALSE
colnames(counts)

# ---- set up deseq2 object ----
#set up experimental design
sample_names <- colnames(counts)

colData <- data.frame(
  sample = sample_names
)

colData <- colData |>
  mutate(
    group = case_when(
      str_detect(sample, "Mock") ~ "Mock",
      TRUE ~ str_extract(sample, "TX|MAD|ROV")
    ),
    time = case_when(
      str_detect(sample, "Mock") ~ "0",
      TRUE ~ str_extract(sample, "12|24|48")
    ),
    condition = ifelse(group == "Mock",
                       "Mock_0",
                       paste0(group, "_", time))
  )

rownames(colData) <- colData$sample
colData$condition <- factor(colData$condition)
colData$condition <- relevel(colData$condition, "Mock_0")



# DESeq model
dds <- DESeqDataSetFromMatrix(
  countData = round(as.matrix(counts)),
  colData   = colData,
  design    = ~ condition
)

dds <- DESeq(dds)

resultsNames(dds)

# ---- Shrink data ----
res <- results(dds)
plotMA(res)


coef_names <- resultsNames(dds)
coef_names <- coef_names[grep("vs_Mock", coef_names)]

res_shrunk_list <- lapply(coef_names, function(coef_name) {
  res <- lfcShrink(
    dds,
    coef = coef_name,
    type = "apeglm"
  )
  res_df <- as.data.frame(res)
  res_df$ensembl_id <- rownames(res_df)
  res_df
})

names(res_shrunk_list) <- coef_names

saveRDS(res_shrunk_list, "data/R2SDHF/Plasmo_res_shrunk_list.rds")