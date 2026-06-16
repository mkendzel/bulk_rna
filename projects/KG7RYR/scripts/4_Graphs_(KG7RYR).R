library(ggplot2)
library(tidyr)
library(dplyr)

# Mouse Ensembl IDs for the genes of interest
genes_of_interest <- c(
  "Nfkb1" = "ENSMUSG00000028163",
  "Irf3"  = "ENSMUSG00000003184"
)

# Define sample groups
groups <- list(
  Control = c("control_liver_1", "control_liver_2", "control_liver_3"),
  CL      = c("cl_liver_4", "cl_liver_5", "cl_liver_6"),
  DOPC    = c("d_liver_7", "d_liver_8", "d_liver_9")
)

# Build a summary data frame with mean counts per group
plot_data <- do.call(rbind, lapply(names(genes_of_interest), function(gene_name) {
  ens_id <- genes_of_interest[[gene_name]]
  do.call(rbind, lapply(names(groups), function(grp) {
    vals <- normalized_counts_liver[ens_id, groups[[grp]]]
    data.frame(
      Gene  = gene_name,
      Group = grp,
      Mean  = mean(vals),
      SE    = sd(vals) / sqrt(length(vals))
    )
  }))
}))

plot_data$Group <- factor(plot_data$Group, levels = c("Control", "CL", "DOPC"))

# Bar plot with SE error bars
ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(y = "Normalized Counts", x = NULL, title = "Liver - Normalized Expression") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

#More genes of interest
genes_of_interest <- c(
  "Ifna1"  = "ENSMUSG00000027397",
  "Ifnb1"  = "ENSMUSG00000048806",
  "Ifng"   = "ENSMUSG00000055170",
  "Il6"    = "ENSMUSG00000025746",
  "Cxcl10" = "ENSMUSG00000034855",
  "Ifit1"  = "ENSMUSG00000034459",
  "Irf7"   = "ENSMUSG00000025785"
)

plot_data <- do.call(rbind, lapply(names(genes_of_interest), function(gene_name) {
  ens_id <- genes_of_interest[[gene_name]]
  do.call(rbind, lapply(names(groups), function(grp) {
    vals <- normalized_counts_liver[ens_id, groups[[grp]]]
    data.frame(
      Gene  = gene_name,
      Group = grp,
      Mean  = mean(vals),
      SE    = sd(vals) / sqrt(length(vals))
    )
  }))
}))

plot_data$Group <- factor(plot_data$Group, levels = c("Control", "CL", "DOPC"))
plot_data$Gene <- factor(plot_data$Gene, levels = names(genes_of_interest))

ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 2) +
  labs(y = "Normalized Counts", x = NULL, title = "Liver - Normalized Expression") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

library(org.Mm.eg.db)

# Just list the gene symbols you want
gene_symbols <- c("Cxcl10", "Ifit1", "Irf7")

# Automatically map symbols to Ensembl IDs
genes_of_interest <- mapIds(org.Mm.eg.db, 
                            keys = gene_symbols, 
                            keytype = "SYMBOL", 
                            column = "ENSEMBL")

# Check for any that didn't map or aren't in the counts matrix
missing <- genes_of_interest[is.na(genes_of_interest) | !genes_of_interest %in% rownames(normalized_counts_liver)]
if (length(missing) > 0) message("Missing genes: ", paste(names(missing), collapse = ", "))

# Filter to only valid genes
genes_of_interest <- genes_of_interest[!is.na(genes_of_interest) & genes_of_interest %in% rownames(normalized_counts_liver)]

# Build plot data
plot_data <- do.call(rbind, lapply(names(genes_of_interest), function(gene_name) {
  ens_id <- genes_of_interest[[gene_name]]
  do.call(rbind, lapply(names(groups), function(grp) {
    vals <- normalized_counts_liver[ens_id, groups[[grp]]]
    data.frame(
      Gene  = gene_name,
      Group = grp,
      Mean  = mean(vals),
      SE    = sd(vals) / sqrt(length(vals))
    )
  }))
}))

plot_data$Group <- factor(plot_data$Group, levels = c("Control", "CL", "DOPC"))
plot_data$Gene <- factor(plot_data$Gene, levels = names(genes_of_interest))

ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1)+
  theme_minimal()+
  labs(y = "Normalized Counts", x = NULL, title = "Liver") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# Define spleen sample groups
groups_spleen <- list(
  Control = c("control_spleen_1", "control_spleen_2", "control_spleen_3"),
  CL      = c("cl_spleen_4", "cl_spleen_5"),
  DOPC    = c("d_spleen_7", "d_spleen_8")
)

# --- Plot 1: Nfkb1 and Irf3 ---
genes_of_interest <- c(
  "Nfkb1" = "ENSMUSG00000028163",
  "Irf3"  = "ENSMUSG00000003184"
)

plot_data <- do.call(rbind, lapply(names(genes_of_interest), function(gene_name) {
  ens_id <- genes_of_interest[[gene_name]]
  do.call(rbind, lapply(names(groups_spleen), function(grp) {
    vals <- normalized_counts_spleen[ens_id, groups_spleen[[grp]]]
    data.frame(
      Gene  = gene_name,
      Group = grp,
      Mean  = mean(vals),
      SE    = sd(vals) / sqrt(length(vals))
    )
  }))
}))

plot_data$Group <- factor(plot_data$Group, levels = c("Control", "CL", "DOPC"))

ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ Gene, scales = "free_y") +
  labs(y = "Normalized Counts", x = NULL, title = "Spleen - Normalized Expression") +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# --- Plot 2: Cxcl10, Ifit1, Irf7 (auto-mapped) ----
library(org.Mm.eg.db)

gene_symbols <- c("Ifng", "Cxcl10", "Ifit1", "Irf7")


genes_of_interest <- mapIds(org.Mm.eg.db, 
                            keys = gene_symbols, 
                            keytype = "SYMBOL", 
                            column = "ENSEMBL")

missing <- genes_of_interest[is.na(genes_of_interest) | !genes_of_interest %in% rownames(normalized_counts_spleen)]
if (length(missing) > 0) message("Missing genes: ", paste(names(missing), collapse = ", "))

genes_of_interest <- genes_of_interest[!is.na(genes_of_interest) & genes_of_interest %in% rownames(normalized_counts_spleen)]

plot_data <- do.call(rbind, lapply(names(genes_of_interest), function(gene_name) {
  ens_id <- genes_of_interest[[gene_name]]
  do.call(rbind, lapply(names(groups_spleen), function(grp) {
    vals <- normalized_counts_spleen[ens_id, groups_spleen[[grp]]]
    data.frame(
      Gene  = gene_name,
      Group = grp,
      Mean  = mean(vals),
      SE    = sd(vals) / sqrt(length(vals))
    )
  }))
}))

plot_data$Group <- factor(plot_data$Group, levels = c("Control", "CL", "DOPC"))
plot_data$Gene <- factor(plot_data$Gene, levels = names(genes_of_interest))

ggplot(plot_data, aes(x = Group, y = Mean, fill = Group)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), width = 0.2) +
  facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
  theme_minimal() +
  labs(y = "Normalized Counts", x = NULL, title = "Spleen") +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

# ---- TF Activity network ----
library(igraph)
library(decoupleR)
library(scales)

# Load the CollecTRI network
tfnet <- get_collectri(organism = "mouse", split_complexes = FALSE)

# Load TF results
tf_results_liver <- readRDS("projects/KG7RYR/data/r_objects/tf_results_liver.rds")

# Pick a contrast
coef_name <- names(tf_results_liver)[2]
tf_res <- tf_results_liver[[coef_name]]

# Top 10 TFs by absolute activity score
tf_top <- tf_res[order(abs(tf_res$score), decreasing = TRUE), ]
tf_top <- head(tf_top, 10)

# Subset network to selected TFs — all targets
net_sub <- tfnet[tfnet$source %in% tf_top$source, ]

# Build graph
g <- graph_from_data_frame(net_sub[, c("source", "target", "mor")], directed = TRUE)

# Color TFs by activity score, targets grey
is_tf <- V(g)$name %in% tf_top$source
tf_score_map <- setNames(tf_top$score, tf_top$source)
tf_scores <- tf_score_map[V(g)$name]
max_abs <- max(abs(tf_scores), na.rm = TRUE)

color_func <- col_numeric(
  palette = c("steelblue", "white", "firebrick"),
  domain = c(-max_abs, max_abs)
)

V(g)$color <- ifelse(is_tf, color_func(tf_scores), "grey80")
V(g)$size <- ifelse(is_tf, 14, 3)
V(g)$label <- ifelse(is_tf, V(g)$name, NA)
V(g)$label.cex <- 0.65
V(g)$frame.color <- NA

# Edge styling
E(g)$width <- 1.5
E(g)$color <- ifelse(E(g)$mor > 0,
                     adjustcolor("firebrick3", alpha.f = 0.3),
                     adjustcolor("steelblue3", alpha.f = 0.3)
)

set.seed(42)
plot(g,
     layout = layout_with_fr(g),
     edge.arrow.mode = 0,
     edge.arrow.size = 0,
     main = paste("TF activity network:", coef_name)
)

# ---- Tables for TF activity ----
tf_res$padj <- p.adjust(tf_res$p_value, method = "BH")
tf_sig <- tf_res[tf_res$padj < 0.05, ]
tf_sig <- tf_sig[order(abs(tf_sig$score), decreasing = TRUE), ]
print(tf_sig, n = nrow(tf_sig))
tf_results_spleen <- readRDS("projects/KG7RYR/data/r_objects/tf_results_spleen.rds")

coef_name <- names(tf_results_spleen)[1]
tf_res <- tf_results_liver[[coef_name]]

tf_res$padj <- p.adjust(tf_res$p_value, method = "BH")
tf_sig <- tf_res[tf_res$padj < 0.05, ]
tf_sig <- tf_sig[order(abs(tf_sig$score), decreasing = TRUE), ]
print(tf_sig, n = nrow(tf_sig))

library(writexl)
write_xlsx(tf_sig, paste0("projects/KG7RYR/data/tf_significant_", coef_name, ".xlsx"))
