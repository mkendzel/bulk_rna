# ==== Library ====
library(fgsea)
library(ggvenn)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(org.Hs.eg.db)  # Mouse: org.Mm.eg.db
library(AnnotationDbi)
library(grid)
library(openxlsx)
library(fgsea)
library(ggrepel)
library(limma)
# ==== Load helper functions ====
invisible(sapply(list.files("R", full.names = TRUE), source))
# ==== Data load ====
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)
mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")
# ---- Brain data load ----
v_brain <- load_checkpoint("voom_brain_d0_d7",
  dir = "projects/mavs_mito/data/r_objects"
)

res_voom_brain_d0_d7_annotated <- load_checkpoint("res_voom_brain_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects"
)

fgsea_wt_vs_mavsko_d7_brain <- load_checkpoint("fgsea_wt_vs_mavsko_d7_brain",
  dir = "projects/mavs_mito/data/r_objects"
)

colData_brain <- load_checkpoint("colData_brain",
  dir = "projects/mavs_mito/data/r_objects"
)

# ======================================
# ==== Brain Graphs ====
# ---- Scatter plot of WT vs MAVS-KO fold changes, colored by significance ----
# +
wt_res <- res_voom_brain_d0_d7_annotated$WT_D7_vs_D0
ko_res <- res_voom_brain_d0_d7_annotated$MAVSKO_D7_vs_D0

# merge WT and MAVS-KO results on ensembl_id, keep relevant columns
merged_fc <- merge(
  wt_res[, c("ensembl_id", "SYMBOL", "logFC", "adj.P.Val", "P.Value")],
  ko_res[, c("ensembl_id", "logFC", "adj.P.Val", "P.Value")],
  by = "ensembl_id",
  suffixes = c("_WT", "_KO")
)

# classify genes by significance (FDR < 0.05) and effect size (|logFC| > 2) in each contrast
sig_wt <- merged_fc$adj.P.Val_WT < 0.05 & abs(merged_fc$logFC_WT) > 2
sig_ko <- merged_fc$adj.P.Val_KO < 0.05 & abs(merged_fc$logFC_KO) > 2

# (optional): fiter
# sig_wt <- merged_fc$P.Value_WT < 0.01 & abs(merged_fc$logFC_WT) > 2
# sig_ko <- merged_fc$P.Value_KO < 0.01 & abs(merged_fc$logFC_KO) > 2

merged_fc$category <- "Not Significant"
merged_fc$category[sig_wt & sig_ko] <- "Both"
merged_fc$category[sig_wt & !sig_ko] <- "WT only"
merged_fc$category[!sig_wt & sig_ko] <- "Mavs-/- only"

merged_fc$category <- factor(
  merged_fc$category,
  levels = c("Not Significant", "Mavs-/- only", "Both", "WT only")
)
merged_fc <- merged_fc[order(merged_fc$category), ]

category_colors <- c(
  "Not Significant" = "grey80",
  "Mavs-/- only" = "firebrick",
  "Both" = "orange",
  "WT only" = "forestgreen"
)

# count genes per category for on-plot labels
counts <- table(merged_fc$category)

x_range <- range(merged_fc$logFC_WT, na.rm = TRUE)
y_range <- range(merged_fc$logFC_KO, na.rm = TRUE)

label_df <- data.frame(
  x = x_range[2] * 0.75,
  y = y_range[2] * 0.9,
  label = paste0("Both\n(", counts["Both"], ")"),
  color = category_colors["Both"]
)

p <- ggplot(merged_fc, aes(x = logFC_WT, y = logFC_KO)) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
  geom_point(aes(fill = category), shape = 21, color = "black", size = 2, stroke = 0.15, alpha = 0.85) +
  geom_text(
    data = label_df,
    aes(x = x, y = y, label = label, color = I(color)),
    fontface = "bold", size = 4, inherit.aes = FALSE
  ) +
  scale_fill_manual(values = category_colors, name = NULL) +
  labs(
    x = "WT Fold Change (log2)",
    y = "Mavs-/- Fold Change (log2)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

p

# directional counts using the significance masks above
# directional breakdown of the plot categories
count_table <- data.frame(
  Category = c(
    "WT only - Up", "WT only - Down",
    "Mavs-/- only - Up", "Mavs-/- only - Down",
    "Both"
  ),
  Genes = c(
    sum(merged_fc$category == "WT only" & merged_fc$logFC_WT > 0),
    sum(merged_fc$category == "WT only" & merged_fc$logFC_WT < 0),
    sum(merged_fc$category == "Mavs-/- only" & merged_fc$logFC_KO > 0),
    sum(merged_fc$category == "Mavs-/- only" & merged_fc$logFC_KO < 0),
    sum(merged_fc$category == "Both")
  )
)
print(count_table)

ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_WT_vs_MAVSKO_FC_scatter.pdf",
  plot     = p,
  width = 4,
  height = 4, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====
# ---- Specific gene scatter plots ----

#+
# merge WT and MAVS-KO results on ensembl_id, keep relevant columns
wt_res <- res_voom_brain_d0_d7_annotated$WT_D7_vs_D0
ko_res <- res_voom_brain_d0_d7_annotated$MAVSKO_D7_vs_D0

merged_fc <- merge(
  wt_res[, c("ensembl_id", "SYMBOL", "logFC", "adj.P.Val", "P.Value")],
  ko_res[, c("ensembl_id", "logFC", "adj.P.Val", "P.Value")],
  by = "ensembl_id",
  suffixes = c("_WT", "_KO")
)

geneset_dir <- "genesets/mavs_mito_genesets"

# read one gene symbol per line, drop blank lines
read_gene_list <- function(path) {
  g <- trimws(readLines(path, warn = FALSE))
  g[nzchar(g)]
}

terminal_effector_genes <- read_gene_list(file.path(geneset_dir, "Goldrath_effector.lists"))
memory_precursor_genes  <- read_gene_list(file.path(geneset_dir, "Goldrathmemory.lists"))

# genes annotated by name on each plot
te_label_genes <- c("Gzmb", "Klrg1", "Prf1", "Ifng")
mp_label_genes <- c("Eomes", "Il15", "Cd44", "Id2")

# case-insensitive membership test against annotated symbols
in_set <- function(symbols, gene_set) toupper(symbols) %in% toupper(gene_set)

# report how many set genes are present in the merged results
report_overlap <- function(df, gene_set, name) {
  n_match <- sum(toupper(gene_set) %in% unique(toupper(df$SYMBOL)))
  message(sprintf("%s: %d of %d genes matched in merged_fc", name, n_match, length(gene_set)))
}

report_overlap(merged_fc, terminal_effector_genes, "Terminal effector")
report_overlap(merged_fc, memory_precursor_genes, "Memory precursor")

plot_signature_scatter <- function(df, sig_genes, label_genes,
                                   highlight_color, x_lab, y_lab) {
  # flag signature membership and draw background points first
  df$highlight <- ifelse(in_set(df$SYMBOL, sig_genes), "sig", "bg")
  df <- df[order(df$highlight), ]

  label_df <- df[in_set(df$SYMBOL, label_genes), ]

  ggplot(df, aes(x = logFC_WT, y = logFC_KO)) +
    geom_hline(yintercept = c(-2, 2), linetype = "solid", color = "grey40", linewidth = 0.4) +
    geom_vline(xintercept = c(-2, 2), linetype = "solid", color = "grey40", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "black", linewidth = 0.5) +
    geom_point(
      data = subset(df, highlight == "bg"),
      fill = "grey80", shape = 21, color = "black", size = 2, stroke = 0.15, alpha = 0.6
    ) +
    geom_point(
      data = subset(df, highlight == "sig"),
      fill = highlight_color, shape = 21, color = "black", size = 2.4, stroke = 0.2, alpha = 0.9
    ) +
    geom_text_repel(
      data = label_df,
      aes(label = SYMBOL),
      size = 3.5, fontface = "italic",
      min.segment.length = 0, segment.color = "grey30", segment.size = 0.3,
      box.padding = 0.5, point.padding = 0.3, max.overlaps = Inf, seed = 1
    ) +
    labs(x = x_lab, y = y_lab) +
    theme_classic(base_size = 14) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      legend.position = "none"
    )
}

p_te <- plot_signature_scatter(
  merged_fc, terminal_effector_genes, te_label_genes,
  highlight_color = "purple3",
  x_lab = "WT Fold Change (log2)", y_lab = "Mavs-/- Fold Change (log2)"
)

p_mp <- plot_signature_scatter(
  merged_fc, memory_precursor_genes, mp_label_genes,
  highlight_color = "deeppink3",
  x_lab = "WT Fold Change (log2)", y_lab = "Mavs-/- Fold Change (log2)"
)

p_te
p_mp
ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_terminal_effector_scatter.pdf",
  plot     = p_te,
  width    = 3, height = 3, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)

ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_memory_precursor_scatter.pdf",
  plot     = p_mp,
  width    = 4, height = 4, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)

# ====

# ---- Effector functions heatmap ----
#+
wt_res <- res_voom_brain_d0_d7_annotated$WT_D7_vs_D0
ko_res <- res_voom_brain_d0_d7_annotated$MAVSKO_D7_vs_D0

# genes for effector functions heatmap, in display order (top to bottom)
effector_genes <- c(
  "Tnfsf10", "Nfatc1", "Tnf", "Il2", "Prf1", "Tnfrsf9", "Fasl", "Ifng", "Gzmb"
)

# pull logFC and adj.P.Val for each gene from both contrasts, using SYMBOL to subset
wt_sub <- wt_res[wt_res$SYMBOL %in% effector_genes, c("SYMBOL", "logFC", "adj.P.Val")]
ko_sub <- ko_res[ko_res$SYMBOL %in% effector_genes, c("SYMBOL", "logFC", "adj.P.Val")]

wt_sub$group <- "WT"
ko_sub$group <- "Mavs-/-"

heatmap_df <- rbind(wt_sub, ko_sub)

# ensure every gene has a row for both groups, filling missing genes with NA
full_grid <- expand.grid(SYMBOL = effector_genes, group = c("WT", "Mavs-/-"))
heatmap_df <- merge(full_grid, heatmap_df, by = c("SYMBOL", "group"), all.x = TRUE)

# set logFC to NA for genes that are not significant,
# so they render as black tiles alongside genes missing from the topTable entirely
heatmap_df$logFC[heatmap_df$adj.P.Val >= 0.05] <- NA

# set factor levels to control row and column order in the plot
heatmap_df$SYMBOL <- factor(heatmap_df$SYMBOL, levels = rev(effector_genes))
heatmap_df$group <- factor(heatmap_df$group, levels = c("WT", "Mavs-/-"))

p <- ggplot(heatmap_df, aes(x = group, y = SYMBOL, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "gray", "yellow", "red"),
    values = scales::rescale(c(-10, -2.50, 0, 2.5, 10)),
    limits = c(-10, 10),
    na.value = "black",
    name = "Fold Change\n(log2)"
  ) +
  labs(title = "Effector\nFunctions", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_effector_functions_heatmap.pdf",
  plot     = p,
  width    = 3, height = 6, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====

# ========================
# ---- Activating receptors heatmap ----
#+
# genes for activating receptors heatmap, in display order (top to bottom)
activating_genes <- c(
  "Itgae", "Il2rg", "Cd27", "Icos", "Cd69",
  "Cd44", "Il2ra", "Klrg1", "Itgax", "Klrk1"
)

# pull logFC and adj.P.Val for each gene from both contrasts, using SYMBOL to subset
wt_sub <- wt_res[wt_res$SYMBOL %in% activating_genes, c("SYMBOL", "logFC", "adj.P.Val")]
ko_sub <- ko_res[ko_res$SYMBOL %in% activating_genes, c("SYMBOL", "logFC", "adj.P.Val")]

wt_sub$group <- "WT"
ko_sub$group <- "Mavs-/-"

heatmap_df <- rbind(wt_sub, ko_sub)

# ensure every gene has a row for both groups, filling missing genes with NA
full_grid <- expand.grid(SYMBOL = activating_genes, group = c("WT", "Mavs-/-"))
heatmap_df <- merge(full_grid, heatmap_df, by = c("SYMBOL", "group"), all.x = TRUE)

# set logFC to NA for genes that are not significant,
# so they render as black tiles alongside genes missing from the topTable entirely
heatmap_df$logFC[heatmap_df$adj.P.Val >= 0.05] <- NA

# set factor levels to control row and column order in the plot
heatmap_df$SYMBOL <- factor(heatmap_df$SYMBOL, levels = rev(activating_genes))
heatmap_df$group <- factor(heatmap_df$group, levels = c("WT", "Mavs-/-"))

p <- ggplot(heatmap_df, aes(x = group, y = SYMBOL, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "gray", "yellow", "red"),
    values = scales::rescale(c(-10, -2.50, 0, 2.5, 10)),
    limits = c(-10, 10),
    na.value = "black",
    name = "Fold Change\n(log2)"
  ) +
  labs(title = "Activating\nReceptors", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_activating_receptors_heatmap.pdf",
  plot     = p,
  width    = 3, height = 6, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====

# ========================
# ---- Cellular stress response heatmap ----
#+
# genes for cellular stress response heatmap, in display order (top to bottom)
stress_genes <- c(
  "Atf4", "Arnt", "Hif1a", "Atf5", "Nt5e", "Atf6",
  "Ddit3", "Eif2ak2", "Hspa1b", "Hspa1a", "Bhlhe40", "Asns"
)

# pull logFC and adj.P.Val for each gene from both contrasts, using SYMBOL to subset
wt_sub <- wt_res[wt_res$SYMBOL %in% stress_genes, c("SYMBOL", "logFC", "adj.P.Val")]
ko_sub <- ko_res[ko_res$SYMBOL %in% stress_genes, c("SYMBOL", "logFC", "adj.P.Val")]

wt_sub$group <- "WT"
ko_sub$group <- "Mavs-/-"

heatmap_df <- rbind(wt_sub, ko_sub)

# ensure every gene has a row for both groups, filling missing genes with NA
full_grid <- expand.grid(SYMBOL = stress_genes, group = c("WT", "Mavs-/-"))
heatmap_df <- merge(full_grid, heatmap_df, by = c("SYMBOL", "group"), all.x = TRUE)

# set logFC to NA for genes that are not significant,
# so they render as black tiles alongside genes missing from the topTable entirely
heatmap_df$logFC[heatmap_df$adj.P.Val >= 0.05] <- NA

# set factor levels to control row and column order in the plot
heatmap_df$SYMBOL <- factor(heatmap_df$SYMBOL, levels = rev(stress_genes))
heatmap_df$group <- factor(heatmap_df$group, levels = c("WT", "Mavs-/-"))

p <- ggplot(heatmap_df, aes(x = group, y = SYMBOL, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "gray", "yellow", "red"),
    values = scales::rescale(c(-10, -2.50, 0, 2.5, 10)),
    limits = c(-10, 10),
    na.value = "black",
    name = "Fold Change\n(log2)"
  ) +
  labs(title = "Cellular Stress\nResponse", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_cellular_stress_response_heatmap.pdf",
  plot     = p,
  width    = 3, height = 6, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====

# ========================
# ---- Mitotic spindle heatmap ----

#+
# genes for mitotic spindle heatmap, in display order (top to bottom)
spindle_genes <- c(
  "Arhgap29", "Tiam1", "Arhgap5", "Kifap3", "Pttg1", "Rapgef5", "Bub1b",
  "Kif22", "Cdk1", "Ccnb1", "Plk1", "Birc5", "Bub1", "Pif1", "Tpx2"
)

# pull logFC and adj.P.Val for each gene from both contrasts, using SYMBOL to subset
wt_sub <- wt_res[wt_res$SYMBOL %in% spindle_genes, c("SYMBOL", "logFC", "adj.P.Val")]
ko_sub <- ko_res[ko_res$SYMBOL %in% spindle_genes, c("SYMBOL", "logFC", "adj.P.Val")]

wt_sub$group <- "WT"
ko_sub$group <- "Mavs-/-"

heatmap_df <- rbind(wt_sub, ko_sub)

# ensure every gene has a row for both groups, filling missing genes with NA
full_grid <- expand.grid(SYMBOL = spindle_genes, group = c("WT", "Mavs-/-"))
heatmap_df <- merge(full_grid, heatmap_df, by = c("SYMBOL", "group"), all.x = TRUE)

# set logFC to NA for genes that are not significant,
# so they render as black tiles alongside genes missing from the topTable entirely
heatmap_df$logFC[heatmap_df$adj.P.Val >= 0.05] <- NA

# set factor levels to control row and column order in the plot
heatmap_df$SYMBOL <- factor(heatmap_df$SYMBOL, levels = rev(spindle_genes))
heatmap_df$group <- factor(heatmap_df$group, levels = c("WT", "Mavs-/-"))

p <- ggplot(heatmap_df, aes(x = group, y = SYMBOL, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradientn(
    colors = c("blue", "lightblue", "gray", "yellow", "red"),
    values = scales::rescale(c(-10, -2.50, 0, 2.5, 10)),
    limits = c(-10, 10),
    na.value = "black",
    name = "Fold Change\n(log2)"
  ) +
  labs(title = "Mitotic\nSpindle", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_text(face = "italic", color = "black"),
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

p
ggsave(
  filename = "projects/mavs_mito/figures/raw/brain_mitotic_spindle_heatmap.pdf",
  plot     = p,
  width    = 3, height = 6, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====
# ---- GSEA curves
# ---- GSEA Curves ----
# ==========================

#+
tt <- res_voom_brain_d0_d7_annotated$WT_vs_MAVSKO_D7
tt <- tt[!is.na(tt$SYMBOL) & tt$SYMBOL != "", ]
tt <- tt[order(tt$AveExpr, decreasing = TRUE), ]
tt <- tt[!duplicated(tt$SYMBOL), ]
ranks <- sort(setNames(tt$t, tt$SYMBOL), decreasing = TRUE)

pathways <- gmtPathways("genesets/mh.all.v2026.1.Mm.symbols.gmt")

# Running enrichment score across the ranked list (GSEA, weighted by |stat|)
running_es <- function(stats, pathway_genes, gsea_param = 1) {
  stats <- sort(stats, decreasing = TRUE)
  n     <- length(stats)
  hits  <- names(stats) %in% pathway_genes
  s     <- abs(stats)^gsea_param
  inc   <- ifelse(hits, s / sum(s[hits]), 0)
  dec   <- ifelse(hits, 0, 1 / (n - sum(hits)))
  es    <- cumsum(inc - dec)
  list(rank = 0:n, es = c(0, es), hit_ranks = which(hits), n = n)
}

bar_cols <- colorRampPalette(c("#3B4CC0", "#F7F7F7", "#B40426"))(256)

# Single pathway panel: curve / gene ticks / gradient bar, stacked
draw_panel <- function(id, title, ranks, pathways, fgsea_res) {
  r   <- running_es(ranks, pathways[[id]])
  i   <- match(id, fgsea_res$pathway)
  lab <- sprintf("NES= %.2f\np-value= %.2g", fgsea_res$NES[i], fgsea_res$pval[i])

  ylim <- range(r$es)
  ypad <- diff(ylim) * 0.1

  layout(matrix(1:3, ncol = 1), heights = c(4, 0.45, 0.7))

  par(mar = c(0, 6, 3.5, 1.5))
  plot(r$rank, r$es, type = "l", lwd = 3, col = "#4477AA",
       xlim = c(0, r$n), ylim = c(ylim[1] - ypad, ylim[2] + ypad),
       axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  abline(h = pretty(ylim), lty = 2, col = "grey75")
  abline(h = 0, col = "black")
  lines(r$rank, r$es, lwd = 3, col = "#4477AA")
  axis(2, las = 1, cex.axis = 1.6)
  mtext("Enrichment Score", side = 2, line = 4, cex = 1.7)
  box()
  title(main = title, font.main = 1, cex.main = 2, line = 1.2)
  text(r$n, ylim[2] + ypad, lab, adj = c(1, 1), cex = 1.6)

  par(mar = c(0, 6, 0, 1.5))
  plot(NA, xlim = c(0, r$n), ylim = c(0, 1), axes = FALSE,
       xlab = "", ylab = "", xaxs = "i")
  segments(r$hit_ranks, 0, r$hit_ranks, 1, lwd = 0.6)

  par(mar = c(4.5, 6, 0, 1.5))
  image(seq(0, r$n, length.out = length(bar_cols)), 1,
        matrix(seq_along(bar_cols), ncol = 1), col = bar_cols,
        axes = FALSE, xlab = "", ylab = "", xaxs = "i")
  box()
  axis(1, cex.axis = 1.6)
  mtext("Rank", side = 1, line = 3, cex = 1.7)
  text(0,   1, "WT",      adj = c(-0.15, 0.5), col = "white", font = 2, cex = 1.6)
  text(r$n, 1, "Mavs-/-", adj = c( 1.15, 0.5), col = "white", font = 2, cex = 1.6)
}
#====

# e2f target
cairo_pdf(
  filename = "projects/mavs_mito/figures/raw/brain_gsea_E2F.pdf",
  width = 6, height = 6
)
draw_panel("HALLMARK_E2F_TARGETS", "E2F Targets",
           ranks, pathways, fgsea_wt_vs_mavsko_d7_brain)
dev.off()


# G2M checkpoint
cairo_pdf(
  filename = "projects/mavs_mito/figures/raw/brain_gsea_G2M.pdf",
  width = 6, height = 6
)
draw_panel("HALLMARK_G2M_CHECKPOINT", "G2M Checkpoint",
           ranks, pathways, fgsea_wt_vs_mavsko_d7_brain)
dev.off()

#oxidative phosphorylation
cairo_pdf(
  filename = "projects/mavs_mito/figures/raw/brain_gsea_OXPHOS.pdf",
  width = 6, height = 6
)
draw_panel("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "Oxidative Phosphorylation",
           ranks, pathways, fgsea_wt_vs_mavsko_d7_brain)
dev.off()

#glycolysis
cairo_pdf(
  filename = "projects/mavs_mito/figures/raw/brain_gsea_GLYCOLYSIS.pdf",
  width = 6, height = 6
)
draw_panel("HALLMARK_GLYCOLYSIS", "Glycolysis",
           ranks, pathways, fgsea_wt_vs_mavsko_d7_brain)
dev.off()

# ==== Spleen ====
# ==== Data load ====
# Human Hallmark gene sets. For mouse, use the matching *.Mm.symbols.gmt file.
gmt_path <- "genesets/h.all.v2025.1.Hs.symbols.gmt"
hallmark_sets <- gmtPathways(gmt_path)
mouse_gene_map <- readRDS("genesets/mouse_gene_map_ensembl_symbol.rds")
v <- load_checkpoint("voom_spleen_d0_d7",
  dir = "projects/mavs_mito/data/r_objects"
)

res_voom_spleen_d0_d7 <- load_checkpoint("res_voom_spleen_d0_d7",
  dir = "projects/mavs_mito/data/r_objects"
)

res_voom_spleen_d0_d7_annotated <- load_checkpoint("res_voom_spleen_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects"
)

fgsea_wt_vs_mavsko_d7 <- load_checkpoint("fgsea_wt_vs_mavsko_d7",
  dir = "projects/mavs_mito/data/r_objects"
)


colData_spleen <- load_checkpoint("colData_spleen",
  dir = "projects/mavs_mito/data/r_objects"
)

# ======================================
# ==== Spleen Graphs ====
# ---- Volcano plot: Mavs-/- vs WT at D0 (naive; WT = reference) ----
# The saved results only include a WT-vs-MAVSKO contrast at D7, so compute the
# equivalent naive (D0) contrast from the saved voom object.
#+
# rebuild the group factor with the reference ordering used when v was fit
colData_spleen$group <- factor(colData_spleen$group,
  levels = c("WT_D0", "MAVSKO_D0", "WT_D7", "MAVSKO_D7")
)
stopifnot(identical(colnames(v), colData_spleen$sample))

design <- model.matrix(~ 0 + group, data = colData_spleen)
colnames(design) <- levels(colData_spleen$group)

contr_d0 <- makeContrasts(WT_vs_MAVSKO_D0 = WT_D0 - MAVSKO_D0, levels = design)
fit_d0 <- eBayes(contrasts.fit(lmFit(v, design), contr_d0))
res_d0 <- topTable(fit_d0, coef = "WT_vs_MAVSKO_D0", number = Inf, sort.by = "none")
res_d0$ensembl_id <- sub("\\..*$", "", rownames(res_d0))

# borrow the ensembl_id -> SYMBOL map already built for this gene universe
sym_map <- res_voom_spleen_d0_d7_annotated$WT_vs_MAVSKO_D7[, c("ensembl_id", "SYMBOL")]
volcano_df <- merge(res_d0, sym_map, by = "ensembl_id")
volcano_df <- volcano_df[!is.na(volcano_df$SYMBOL) & volcano_df$SYMBOL != "", ]

# contrast is WT - MAVSKO, so flip the sign to make positive = up in Mavs-/-
volcano_df$logFC_vs_WT <- -volcano_df$logFC
# use FDR on the y-axis to match the FDR-based significance categories below
volcano_df$neglog10FDR <- -log10(volcano_df$adj.P.Val)

# significance thresholds
fc_cut <- 1.5        # |log2 fold change| cutoff
fdr_cut <- 0.05    # adjusted p-value cutoff

volcano_df$category <- "Not Significant"
volcano_df$category[volcano_df$adj.P.Val < fdr_cut & volcano_df$logFC_vs_WT >  fc_cut] <- "Up in Mavs-/-"
volcano_df$category[volcano_df$adj.P.Val < fdr_cut & volcano_df$logFC_vs_WT < -fc_cut] <- "Up in WT"

volcano_df$category <- factor(
  volcano_df$category,
  levels = c("Not Significant", "Up in WT", "Up in Mavs-/-")
)
volcano_df <- volcano_df[order(volcano_df$category), ]

volcano_colors <- c(
  "Not Significant" = "grey80",
  "Up in WT" = "forestgreen",
  "Up in Mavs-/-" = "firebrick"
)

# report how many genes pass the significance thresholds (FDR < fdr_cut & |logFC| > fc_cut)
sig_counts <- table(volcano_df$category)
message(sprintf(
  "Significant genes (FDR < %g, |log2FC| > %g): %d total (%d Up in WT, %d Up in Mavs-/-)",
  fdr_cut, fc_cut,
  sig_counts["Up in WT"] + sig_counts["Up in Mavs-/-"],
  sig_counts["Up in WT"], sig_counts["Up in Mavs-/-"]
))

# label the top genes by significance in each direction
n_label <- 15
sig_df <- volcano_df[volcano_df$category != "Not Significant", ]
label_df <- head(sig_df[order(sig_df$adj.P.Val), ], n_label)

p <- ggplot(volcano_df, aes(x = logFC_vs_WT, y = neglog10FDR)) +
  geom_vline(xintercept = c(-fc_cut, fc_cut), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_point(aes(fill = category), shape = 21, color = "black", size = 2, stroke = 0.15, alpha = 0.85) +
  geom_text_repel(
    data = label_df,
    aes(label = SYMBOL),
    size = 3, fontface = "italic",
    min.segment.length = 0, segment.color = "grey30", segment.size = 0.3,
    box.padding = 0.4, point.padding = 0.3, max.overlaps = Inf, seed = 1
  ) +
  scale_fill_manual(values = volcano_colors, name = NULL) +
  labs(
    x = "Fold Change over WT (log2)",
    y = expression(-log[10]~"FDR"),
    title = "Spleen: WT vs Mavs-/- at D0 (naive)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none"
  )

p
ggsave(
  filename = "projects/mavs_mito/figures/raw/spleen_WT_vs_MAVSKO_D0_volcano_nolegend.pdf",
  plot     = p,
  width    = 5, height = 5, units = "in",
  dpi      = 300,
  device   = cairo_pdf
)
# ====