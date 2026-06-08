# ---- Libraries ----
library(decoupleR)
library(OmnipathR)
library(ggplot2)
# ---- Data Load ----
# Get mouse DoRothEA network or Collec
tfnet <- get_dorothea(organism = "mouse", levels = c("A", "B", "C"))
tfnet <- get_collectri(organism = "mouse", split_complexes = FALSE)

# Full interactive model results
res_shrunk_list <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_qcmin10_annotated.rds")
dds <- readRDS("data/KG7RYR/r_objects/plasmo_dds_raw.rds")

coef_names <- DESeq2::resultsNames(dds)[-1]
# Extract unshrunken results which contain the Wald statistic and p-values
res_wald_list <- setNames(
  lapply(coef_names, function(coef_name) {
    DESeq2::results(dds, name = coef_name)
  }),
  coef_names
)

# liver Results
res_shrunk_liver <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_liver_qcmin10_annotated.rds")
res_wald_liver <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_liver_qcmin10_annotated.rds")

# Spleen Results
res_wald_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resWald_spleen_qcmin10_annotated.rds")
res_shrunk_spleen <- readRDS("data/KG7RYR/r_objects/plasmo_resShrink_spleen_qcmin10_annotated.rds")

# =========================================
# ---- DoRothEA: Transcription factors ----
# =========================================

# ---- Full Interactive model ----
# Run decoupleR using gene symbols from the shrunk results
tf_results_list <- setNames(
  lapply(coef_names, function(coef_name) {
    res <- res_wald_list[[coef_name]]
    shrunk <- res_shrunk_list[[coef_name]]
    
    stats <- res$stat
    symbols <- shrunk$symbol
    names(stats) <- symbols
    
    valid <- !is.na(stats) & !is.na(symbols) & symbols != ""
    stats <- stats[valid]
    
    stats <- tapply(stats, names(stats), function(x) x[which.max(abs(x))])
    stats <- as.numeric(stats)
    names(stats) <- names(tapply(res$stat[valid], symbols[valid], function(x) x[which.max(abs(x))]))
    
    # Genes as rows, contrast as column
    stat_mat <- matrix(stats, ncol = 1, dimnames = list(names(stats), coef_name))
    
    run_ulm(
      mat = stat_mat,
      network = tfnet,
      .source = source,
      .target = target,
      minsize = 5
    )
  }),
  coef_names
)

#Save/load
saveRDS(tf_results_list, "data/KG7RYR/r_objects/tf_results_list.rds")
tf_results_list <- readRDS("data/KG7RYR/r_objects/tf_results_list.rds")

tf_results_all <- do.call(rbind, tf_results_list)
tf_consensus <- tf_results_all[tf_results_all$statistic == "consensus", ]

#Graph

lapply(coef_names, function(coef_name) {
  tf <- tf_results_list[[coef_name]]
  tf <- tf[order(abs(tf$score), decreasing = TRUE), ]
  tf <- head(tf, 25)
  tf$source <- reorder(tf$source, tf$score)
  
  ggplot(tf, aes(x = score, y = source, fill = score > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
    labs(title = coef_name, x = "TF activity (ULM score)", y = NULL) +
    theme_minimal()
})

# ---- DoRothEA - liver ----
coef_names_liver <- names(res_wald_liver)

tf_results_liver <- setNames(
  lapply(coef_names_liver, function(coef_name) {
    res <- res_wald_liver[[coef_name]]
    
    stats <- res$stat
    symbols <- res$symbol
    names(stats) <- symbols
    
    valid <- !is.na(stats) & !is.na(symbols) & symbols != ""
    stats <- stats[valid]
    
    stats <- tapply(stats, names(stats), function(x) x[which.max(abs(x))])
    stats <- as.numeric(stats)
    names(stats) <- names(tapply(res$stat[valid], res$symbol[valid], function(x) x[which.max(abs(x))]))
    
    stat_mat <- matrix(stats, ncol = 1, dimnames = list(names(stats), coef_name))
    
    run_ulm(
      mat = stat_mat,
      network = tfnet,
      .source = source,
      .target = target,
      minsize = 5
    )
  }),
  coef_names_liver
)
#Save
saveRDS(tf_results_liver, "data/KG7RYR/r_objects/tf_results_liver.rds")
tf_results_liver <- readRDS("data/KG7RYR/r_objects/tf_results_liver.rds")


lapply(names(tf_results_liver), function(coef_name) {
  tf <- tf_results_liver[[coef_name]]
  tf <- tf[order(abs(tf$score), decreasing = TRUE), ]
  tf <- head(tf, 25)
  tf$source <- reorder(tf$source, tf$score)
  
  ggplot(tf, aes(x = score, y = source, fill = score > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
    labs(title = paste("liver -", coef_name), x = "TF activity (ULM score)", y = NULL) +
    theme_minimal()
})


tf <- tf_results_liver[[3]]
tf <- tf[order(abs(tf$score), decreasing = TRUE), ]
tf <- head(tf, 25)
print(tf[, c("source", "score", "p_value")], n = 25)
# ---- TF Activity - Spleen ----
coef_names_spleen <- names(res_wald_spleen)

tf_results_spleen <- setNames(
  lapply(coef_names_spleen, function(coef_name) {
    res <- res_wald_spleen[[coef_name]]
    
    stats <- res$stat
    symbols <- res$symbol
    names(stats) <- symbols
    
    valid <- !is.na(stats) & !is.na(symbols) & symbols != ""
    stats <- stats[valid]
    
    stats <- tapply(stats, names(stats), function(x) x[which.max(abs(x))])
    stats <- as.numeric(stats)
    names(stats) <- names(tapply(res$stat[valid], res$symbol[valid], function(x) x[which.max(abs(x))]))
    
    stat_mat <- matrix(stats, ncol = 1, dimnames = list(names(stats), coef_name))
    
    run_ulm(
      mat = stat_mat,
      network = tfnet,
      .source = source,
      .target = target,
      minsize = 5
    )
  }),
  coef_names_spleen
)

#
saveRDS(tf_results_spleen, "data/KG7RYR/r_objects/tf_results_spleen.rds")
tf_results_spleen <- readRDS("data/KG7RYR/r_objects/tf_results_spleen.rds")

lapply(names(tf_results_spleen), function(coef_name) {
  tf <- tf_results_spleen[[coef_name]]
  tf <- tf[order(abs(tf$score), decreasing = TRUE), ]
  tf <- head(tf, 25)
  tf$source <- reorder(tf$source, tf$score)
  
  ggplot(tf, aes(x = score, y = source, fill = score > 0)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"), guide = "none") +
    labs(title = paste("Spleen -", coef_name), x = "TF activity (ULM score)", y = NULL) +
    theme_minimal()
})

# ---- Save all TF results as excel ----
library(openxlsx)

# Define tissues and their corresponding list objects
tissues <- list(
  liver = tf_results_liver,
  spleen = tf_results_spleen
)

# Define the comparisons and their clean names
comparisons <- c(
  "treatment_cl_vs_control" = "CL",
  "treatment_d_vs_control" = "DOPC"
)

# Adjust p-values, filter for significance, and write to Excel
for (tissue in names(tissues)) {
  for (comp in names(comparisons)) {
    df <- tissues[[tissue]][[comp]]
    df$padj <- p.adjust(df$p_value, method = "BH")
    sig <- df[df$padj < 0.05 & !is.na(df$padj), ]
    filename <- file.path("data", "KG7RYR", paste0(tissue, "_", comparisons[comp], ".xlsx"))
    write.xlsx(sig, file = filename)
  }
}
