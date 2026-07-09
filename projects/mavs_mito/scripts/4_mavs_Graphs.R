# ---- Library ----
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
# ==== Load helper functions ====
invisible(sapply(list.files("R", full.names = TRUE), source))
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
# ---- Lolipop plots of GSEA results ----
 
# ---- Average expresion bar graphs ----
# =====================================

# Pull Mavs stats from the three relevant contrasts
get_mavs_stat <- function(contrast_name) {
  df <- res_voom_spleen_d0_d7_annotated[[contrast_name]]
  df[!is.na(df$SYMBOL) & df$SYMBOL == "Mavs", c("SYMBOL", "logFC", "adj.P.Val")]
}

mavs_ko_day    <- get_mavs_stat("MAVSKO_D7_vs_D0")

sig_label <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "NS"
}

mavs_ensembl <- mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL == "Mavs"]
mavs_row     <- grep(paste(mavs_ensembl, collapse = "|"), rownames(v$E), value = TRUE)

df <- data.frame(
  log2cpm = as.numeric(v$E[mavs_row, ]),
  group   = factor(colData_spleen$group,
                   levels = c("MAVSKO_D0", "MAVSKO_D7"))
)
df <- df[!is.na(df$group), ]

group_stats <- do.call(rbind, lapply(split(df$log2cpm, df$group), function(x) {
  data.frame(mean = mean(x), sd = sd(x))
}))
group_stats$group <- factor(rownames(group_stats),
                            levels = c("MAVSKO_D0", "MAVSKO_D7"))

ymax <- max(group_stats$mean + group_stats$sd)
y1 <- ymax * 1.05   # MAVSKO D0 vs D7 bracket

mavs_plot <- ggplot(group_stats, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, linewidth = 0.8) +
  geom_jitter(data = df, aes(x = group, y = log2cpm, fill = group),
              width = 0.1, size = 2, shape = 21, inherit.aes = FALSE) +
  # MAVSKO D0 vs D7
  annotate("segment", x = 1, xend = 2, y = y1, yend = y1) +
  annotate("text", x = 1.5, y = y1 * 1.05,
           label = sig_label(mavs_ko_day$adj.P.Val), size = 5) +
  scale_fill_manual(values = c("MAVSKO_D0" = "grey70", "MAVSKO_D7" = "grey30")) +
  scale_x_discrete(labels = c("MAVS -/- D0", "MAVS -/- D7")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Normalized Read Count (log2)", title = expression(italic("Mavs"))) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# save
ggsave(filename = "projects/mavs_mito/figures/raw/spleen_mavs_barplot.pdf",
       plot     = mavs_plot,
       width    = 3, height = 4, units = "in",
       dpi      = 300,
       device   = cairo_pdf)
# ========================
# ---- Specific genes bar graphs----
# ========================

genes_of_interest <- c("Rigi", "Dhx58", "Ifih1", "Traf6", "Tbk1",
                        "Irf3", "Irf7", "Ifnb1", "Ifna2")

# Display names for facet strips (internal symbol -> paper name)
display_names <- c(Rigi = "Ddx58", Dhx58 = "Dhx58", Ifih1 = "Ifih1",
                   Traf6 = "Traf6", Tbk1 = "Tbk1", Irf3 = "Irf3",
                   Irf7 = "Irf7", Ifnb1 = "Ifnb1", Ifna2 = "Ifna2")

# Resolve each symbol to a single v$E row, highest mean expression if multiple IDs map
get_expr_row <- function(sym) {
  ens_ids <- mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL == sym]
  if (length(ens_ids) == 0) return(NULL)
  rows <- grep(paste(ens_ids, collapse = "|"), rownames(v$E), value = TRUE)
  if (length(rows) == 0) return(NULL)
  if (length(rows) > 1) rows <- rows[which.max(rowMeans(v$E[rows, , drop = FALSE]))]
  rows
}

df_long <- do.call(rbind, lapply(genes_of_interest, function(sym) {
  row <- get_expr_row(sym)
  if (is.null(row)) return(NULL)
  data.frame(
    SYMBOL  = sym,
    log2cpm = as.numeric(v$E[row, ]),
    group   = factor(colData_spleen$group,
                     levels = c("MAVSKO_D0", "MAVSKO_D7"))
  )
}))
df_long <- df_long[!is.na(df_long$group), ]

# Keep ALL genes as x-axis slots; undetected ones render as empty positions
df_long$SYMBOL <- factor(df_long$SYMBOL, levels = genes_of_interest)

# Group means and SDs (only for detected genes)
group_stats <- do.call(rbind, lapply(split(df_long, list(df_long$SYMBOL, df_long$group), drop = TRUE), function(x) {
  data.frame(SYMBOL = x$SYMBOL[1], group = x$group[1],
             mean = mean(x$log2cpm), sd = sd(x$log2cpm))
}))
group_stats$SYMBOL <- factor(group_stats$SYMBOL, levels = genes_of_interest)
group_stats$group  <- factor(group_stats$group,
                             levels = c("MAVSKO_D0", "MAVSKO_D7"))

detected_genes <- levels(droplevels(group_stats$SYMBOL))

# Pull adj.P.Val for each gene from the MAVSKO day contrast
get_pval <- function(contrast, sym) {
  df   <- res_voom_spleen_d0_d7_annotated[[contrast]]
  vals <- df$adj.P.Val[!is.na(df$SYMBOL) & df$SYMBOL == sym]
  if (length(vals) == 0) return(NA)
  vals[1]
}

sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

# One D0-vs-D7 bracket + star per detected gene, centered over its two dodged bars
dodge_w    <- 0.7
bar_offset <- dodge_w / 4   # half-distance between the two dodged bars

sig_df <- do.call(rbind, lapply(detected_genes, function(sym) {
  sub  <- group_stats[group_stats$SYMBOL == sym, ]
  ymax <- max(sub$mean + sub$sd, na.rm = TRUE)
  xpos <- match(sym, genes_of_interest)
  data.frame(
    SYMBOL = sym,
    label  = sig_label(get_pval("MAVSKO_D7_vs_D0", sym)),
    x      = xpos - bar_offset,
    xend   = xpos + bar_offset,
    xmid   = xpos,
    y      = ymax * 1.08
  )
}))
sig_df$SYMBOL <- factor(sig_df$SYMBOL, levels = genes_of_interest)

group_colors <- c("MAVSKO_D0" = "grey75", "MAVSKO_D7" = "grey35")
pd <- position_dodge(width = dodge_w)

ggplot(group_stats, aes(x = SYMBOL, y = mean, fill = group)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.4, position = pd) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.18, linewidth = 0.6, position = pd) +
  geom_point(data = df_long, aes(x = SYMBOL, y = log2cpm, fill = group),
             position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_w),
             size = 1.6, shape = 21, stroke = 0.3, color = "black",
             alpha = 0.9, inherit.aes = FALSE) +
  geom_segment(data = sig_df, aes(x = x, xend = xend, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = xmid, y = y, label = label),
            inherit.aes = FALSE, size = 4, vjust = -0.2) +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = display_names, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = expression("Normalized Read Count (log"[2]*")")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x     = element_text(face = "italic", size = 12, angle = 45, hjust = 1),
    axis.text.y     = element_text(size = 10)
  )

ggsave(filename = "projects/mavs_mito/figures/raw/spleen_genes_barplot.pdf",
       width    = 6, height = 4, units = "in",
       dpi      = 300,
       device   = cairo_pdf)
# ==========================
# ---- GSEA Curves ----
# ==========================


tt <- res_voom_spleen_d0_d7_annotated$WT_vs_MAVSKO_D7
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

# Call once per pathway, one figure at a time
draw_panel("HALLMARK_E2F_TARGETS", "E2F Targets",
           ranks, pathways, fgsea_wt_vs_mavsko_d7)


# ---- Scatter plot of WT vs MAVS-KO fold changes, colored by significance ----
# +


# ====

# ========================
# ---- Cellular stress response heatmap ----
#+
wt_res <- res_voom_spleen_d0_d7_annotated$WT_D7_vs_D0
ko_res <- res_voom_spleen_d0_d7_annotated$MAVSKO_D7_vs_D0

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

# ###########################################################################
# ==== BRAIN GRAPHS ====
# ###########################################################################
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
# ---- Brain: Average expresion bar graphs ----
# =====================================

# Pull Mavs stats from the three relevant contrasts
get_mavs_stat <- function(contrast_name) {
  df <- res_voom_brain_d0_d7_annotated[[contrast_name]]
  df[!is.na(df$SYMBOL) & df$SYMBOL == "Mavs", c("SYMBOL", "logFC", "adj.P.Val")]
}

mavs_wt_day    <- get_mavs_stat("WT_D7_vs_D0")
mavs_ko_day    <- get_mavs_stat("MAVSKO_D7_vs_D0")
mavs_genotype  <- get_mavs_stat("WT_vs_MAVSKO_D7")

sig_label <- function(p) {
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "NS"
}

mavs_ensembl <- mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL == "Mavs"]
mavs_row     <- grep(paste(mavs_ensembl, collapse = "|"), rownames(v_brain$E), value = TRUE)

df <- data.frame(
  log2cpm = as.numeric(v_brain$E[mavs_row, ]),
  group   = factor(colData_brain$group,
                   levels = c("WT_D0", "WT_D7", "MAVSKO_D0", "MAVSKO_D7"))
)

group_stats <- do.call(rbind, lapply(split(df$log2cpm, df$group), function(x) {
  data.frame(mean = mean(x), sd = sd(x))
}))
group_stats$group <- factor(rownames(group_stats),
                            levels = c("WT_D0", "WT_D7", "MAVSKO_D0", "MAVSKO_D7"))

ymax <- max(group_stats$mean + group_stats$sd)
y1 <- ymax * 1.05   # WT D0 vs D7 bracket
y2 <- ymax * 1.15   # MAVSKO D0 vs D7 bracket
y3 <- ymax * 1.25   # WT D7 vs MAVSKO D7 bracket

ggplot(group_stats, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.2, linewidth = 0.8) +
  geom_jitter(data = df, aes(x = group, y = log2cpm, fill = group),
              width = 0.1, size = 2, shape = 21, inherit.aes = FALSE) +
  # WT D0 vs D7
  annotate("segment", x = 1, xend = 2, y = y1, yend = y1) +
  annotate("text", x = 1.5, y = y1 * 1.01,
           label = sig_label(mavs_wt_day$adj.P.Val), size = 5) +
  # MAVSKO D0 vs D7
  annotate("segment", x = 3, xend = 4, y = y2, yend = y2) +
  annotate("text", x = 3.5, y = y2 * 1.01,
           label = sig_label(mavs_ko_day$adj.P.Val), size = 5) +
  # WT D7 vs MAVSKO D7
  annotate("segment", x = 2, xend = 4, y = y3, yend = y3) +
  annotate("text", x = 3, y = y3 * 1.01,
           label = sig_label(mavs_genotype$adj.P.Val), size = 5) +
  scale_fill_manual(values = c("WT_D0" = "grey70", "WT_D7" = "grey30",
                               "MAVSKO_D0" = "#6699CC", "MAVSKO_D7" = "#2255AA")) +
  scale_x_discrete(labels = c("WT D0", "WT D7", "MAVSKO D0", "MAVSKO D7")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Normalized Read Count (log2)", title = expression(italic("Mavs"))) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5))

# ========================
# ---- Brain: Specific genes bar graphs----
# ========================

genes_of_interest <- c("Rigi", "Dhx58", "Ifih1", "Traf6", "Tbk1",
                        "Irf3", "Irf7", "Ifnb1", "Ifna2")

# Display names for facet strips (internal symbol -> paper name)
display_names <- c(Rigi = "Ddx58", Dhx58 = "Dhx58", Ifih1 = "Ifih1",
                   Traf6 = "Traf6", Tbk1 = "Tbk1", Irf3 = "Irf3",
                   Irf7 = "Irf7", Ifnb1 = "Ifnb1", Ifna2 = "Ifna2")

# Resolve each symbol to a single v_brain$E row, highest mean expression if multiple IDs map
get_expr_row <- function(sym) {
  ens_ids <- mouse_gene_map$ENSEMBL[mouse_gene_map$SYMBOL == sym]
  if (length(ens_ids) == 0) return(NULL)
  rows <- grep(paste(ens_ids, collapse = "|"), rownames(v_brain$E), value = TRUE)
  if (length(rows) == 0) return(NULL)
  if (length(rows) > 1) rows <- rows[which.max(rowMeans(v_brain$E[rows, , drop = FALSE]))]
  rows
}

df_long <- do.call(rbind, lapply(genes_of_interest, function(sym) {
  row <- get_expr_row(sym)
  if (is.null(row)) return(NULL)
  data.frame(
    SYMBOL  = sym,
    log2cpm = as.numeric(v_brain$E[row, ]),
    group   = factor(colData_brain$group,
                     levels = c("WT_D0", "WT_D7", "MAVSKO_D0", "MAVSKO_D7"))
  )
}))

# Keep ALL genes as factor levels so undetected ones render as empty facets
df_long$SYMBOL <- factor(df_long$SYMBOL, levels = genes_of_interest)

# Group means and SDs (only for detected genes)
group_stats <- do.call(rbind, lapply(split(df_long, list(df_long$SYMBOL, df_long$group), drop = TRUE), function(x) {
  data.frame(SYMBOL = x$SYMBOL[1], group = x$group[1],
             mean = mean(x$log2cpm), sd = sd(x$log2cpm))
}))
group_stats$SYMBOL <- factor(group_stats$SYMBOL, levels = genes_of_interest)
group_stats$group  <- factor(group_stats$group,
                             levels = c("WT_D0", "WT_D7", "MAVSKO_D0", "MAVSKO_D7"))

detected_genes <- levels(droplevels(group_stats$SYMBOL))

# Pull adj.P.Val for each gene from the three contrasts
get_pval <- function(contrast, sym) {
  df   <- res_voom_brain_d0_d7_annotated[[contrast]]
  vals <- df$adj.P.Val[!is.na(df$SYMBOL) & df$SYMBOL == sym]
  if (length(vals) == 0) return(NA)
  vals[1]
}

sig_label <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) "***"
  else if (p < 0.01) "**"
  else if (p < 0.05) "*"
  else "ns"
}

# Brackets only for detected genes
sig_df <- do.call(rbind, lapply(detected_genes, function(sym) {
  data.frame(
    SYMBOL      = sym,
    wt_day      = sig_label(get_pval("WT_D7_vs_D0",     sym)),
    ko_day      = sig_label(get_pval("MAVSKO_D7_vs_D0", sym)),
    genotype_d7 = sig_label(get_pval("WT_vs_MAVSKO_D7", sym))
  )
}))

ypos <- do.call(rbind, lapply(detected_genes, function(sym) {
  sub  <- group_stats[group_stats$SYMBOL == sym, ]
  ymax <- max(sub$mean + sub$sd, na.rm = TRUE)
  data.frame(SYMBOL = sym, y1 = ymax * 1.10, y2 = ymax * 1.22, y3 = ymax * 1.34)
}))

sig_df <- merge(sig_df, ypos, by = "SYMBOL")
sig_df$SYMBOL <- factor(sig_df$SYMBOL, levels = genes_of_interest)

group_colors <- c("WT_D0" = "grey75", "WT_D7" = "grey35",
                  "MAVSKO_D0" = "#9DBEDC", "MAVSKO_D7" = "#2F5C9E")

ggplot(group_stats, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.4) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.18, linewidth = 0.6) +
  geom_jitter(data = df_long, aes(x = group, y = log2cpm, fill = group),
              width = 0.12, size = 1.6, shape = 21, stroke = 0.3,
              color = "black", alpha = 0.9, inherit.aes = FALSE) +
  geom_segment(data = sig_df, aes(x = 1, xend = 2, y = y1, yend = y1),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = 1.5, y = y1, label = wt_day),
            inherit.aes = FALSE, size = 4, vjust = -0.2) +
  geom_segment(data = sig_df, aes(x = 3, xend = 4, y = y2, yend = y2),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = 3.5, y = y2, label = ko_day),
            inherit.aes = FALSE, size = 4, vjust = -0.2) +
  geom_segment(data = sig_df, aes(x = 2, xend = 4, y = y3, yend = y3),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = sig_df, aes(x = 3, y = y3, label = genotype_d7),
            inherit.aes = FALSE, size = 4, vjust = -0.2) +
  scale_fill_manual(values = group_colors) +
  scale_x_discrete(labels = c("WT\nD0", "WT\nD7", "KO\nD0", "KO\nD7"),
                   drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  facet_wrap(~ SYMBOL, scales = "free_y", nrow = 2, drop = FALSE,
             labeller = as_labeller(display_names)) +
  labs(x = NULL, y = expression("Normalized Read Count (log"[2]*")")) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    strip.text      = element_text(face = "italic", size = 14),
    strip.background = element_blank(),
    axis.text.x     = element_text(size = 10),
    axis.text.y     = element_text(size = 10),
    panel.spacing   = unit(1, "lines")
  )

# ==========================
# ---- Brain: GSEA Curves ----
# ==========================
# Reuses running_es(), draw_panel(), bar_cols, and pathways from the spleen section.

tt <- res_voom_brain_d0_d7_annotated$WT_vs_MAVSKO_D7
tt <- tt[!is.na(tt$SYMBOL) & tt$SYMBOL != "", ]
tt <- tt[order(tt$AveExpr, decreasing = TRUE), ]
tt <- tt[!duplicated(tt$SYMBOL), ]
ranks <- sort(setNames(tt$t, tt$SYMBOL), decreasing = TRUE)

# Call once per pathway, one figure at a time
draw_panel("HALLMARK_E2F_TARGETS", "E2F Targets",
           ranks, pathways, fgsea_wt_vs_mavsko_d7_brain)

# ========================
# ---- Brain: Cellular stress response heatmap ----
# ========================
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

# ====
# ---- Template plots ----
# ========================

# ---- PCA ----
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_data$treatment_time <- as.character(pca_data$condition)

# EDIT: replace with your own condition/treatment names and colors.
# TreatmentB is omitted since it was dropped during QC in script 1.
base_cols <- c(
  Control     = "#D62728",
  TreatmentA  = "blue",
  TreatmentB  = "darkgreen",
  TreatmentC  = "purple"
)

make_shades <- function(hex, n = 3) {
  colorRampPalette(c("#FFFFFF", hex))(n + 1)[-1]
}

a_shades <- make_shades(base_cols[["TreatmentA"]], 3)
c_shades <- make_shades(base_cols[["TreatmentC"]], 3)

col_map <- c(
  Control_0    = base_cols[["Control"]],
  TreatmentA_1 = a_shades[1], TreatmentA_2 = a_shades[2], TreatmentA_3 = a_shades[3],
  TreatmentC_1 = c_shades[1], TreatmentC_2 = c_shades[2], TreatmentC_3 = c_shades[3]
)

pca_data$treatment_time <- factor(
  pca_data$treatment_time,
  levels = c("Control_0","TreatmentA_1","TreatmentA_2","TreatmentA_3",
             "TreatmentC_1","TreatmentC_2","TreatmentC_3")
)


centroids <- pca_data %>%
  group_by(treatment_time) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    .groups = "drop"
  )

p <- ggplot(pca_data, aes(PC1, PC2, color = treatment_time)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_point(
    data = centroids,
    aes(PC1, PC2, fill = treatment_time),
    shape = 24,           # filled triangle
    color = "black",
    size = 3.5,
    stroke = 1.2
  ) +
  scale_color_manual(values = col_map, drop = FALSE) +
  scale_fill_manual(values = col_map, drop = FALSE) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_classic()
p

ggsave(
  filename = "projects/PROJECT_NAME/figures/PCA_treatment_time.png",
  plot = p,
  width = 8,
  height = 6
)

# ---- Volcano Plots ----
# EDIT: match the significance cutoffs used elsewhere
padj_cutoff <- 0.05
lfc_cutoff  <- 1.5

outdir <- "projects/PROJECT_NAME/figures/volcano"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (coef_name in names(res_shrunk_list)) {
  
  volcano_df <- as.data.frame(res_shrunk_list[[coef_name]]) |>
    tibble::rownames_to_column("ensembl_id") |>
    tibble::as_tibble() |>
    dplyr::mutate(
      neglog10_padj = -log10(padj),
      direction = dplyr::case_when(
        padj < padj_cutoff & log2FoldChange >  lfc_cutoff ~ "Up",
        padj < padj_cutoff & log2FoldChange < -lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  n_up   <- sum(volcano_df$direction == "Up", na.rm = TRUE)
  n_down <- sum(volcano_df$direction == "Down", na.rm = TRUE)
  
  x_max <- max(volcano_df$log2FoldChange, na.rm = TRUE)
  x_min <- min(volcano_df$log2FoldChange, na.rm = TRUE)
  y_max <- max(volcano_df$neglog10_padj[is.finite(volcano_df$neglog10_padj)], na.rm = TRUE)
  
  p <- ggplot2::ggplot(volcano_df, ggplot2::aes(x = log2FoldChange, y = neglog10_padj)) +
    ggplot2::geom_point(ggplot2::aes(color = direction), alpha = 0.6, size = 1) +
    ggplot2::geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed") +
    ggplot2::annotate(
      "text", x = x_max, y = y_max,
      label = paste0("Up: ", n_up),
      hjust = 1, vjust = 1, color = "red"
    ) +
    ggplot2::annotate(
      "text", x = x_min, y = y_max,
      label = paste0("Down: ", n_down),
      hjust = 0, vjust = 1, color = "blue"
    ) +
    ggplot2::scale_color_manual(values = c("Up" = "red", "Down" = "blue", "NS" = "grey70")) +
    ggplot2::theme_classic() +
    ggplot2::labs(
      title = coef_name,
      x = "log2 Fold Change (apeglm-shrunk)",
      y = "-log10 adjusted p-value"
    )
  
  comp_label <- sub("^condition_", "", coef_name)
  outfile <- file.path(outdir, paste0("volcano_", comp_label, ".png"))
  
  ggplot2::ggsave(outfile, p, width = 7, height = 6, dpi = 300)
}

# ---- Heatmaps highly variable ----

vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)

# Ensembl -> SYMBOL map via org.Hs.eg.db
ens_ids <- rownames(mat)

map_df <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys     = ens_ids,
  keytype  = "ENSEMBL",
  columns  = c("SYMBOL")
) |>
  tibble::as_tibble() |>
  dplyr::filter(!is.na(SYMBOL) & SYMBOL != "") |>
  dplyr::distinct(ENSEMBL, .keep_all = TRUE) |>
  dplyr::rename(ensembl_id = ENSEMBL, gene_name = SYMBOL)

idx <- match(ens_ids, map_df$ensembl_id)
mapped_names <- map_df$gene_name[idx]

# Top 50 most variable genes
gene_vars <- apply(mat, 1, var, na.rm = TRUE)
top_ens <- names(sort(gene_vars, decreasing = TRUE))[1:50]

mat_subset <- mat[top_ens, , drop = FALSE]

# Replace Ensembl rownames with gene_name (fallback to Ensembl)
new_names <- mapped_names[match(top_ens, ens_ids)]
new_names[is.na(new_names) | new_names == ""] <- top_ens[is.na(new_names) | new_names == ""]
rownames(mat_subset) <- new_names

# Z-score per gene
mat_scaled <- t(scale(t(mat_subset)))

# Column annotations from colData(dds)

cd <- as.data.frame(colData(dds))

if (!"sample" %in% colnames(cd)) {
  cd$sample <- rownames(cd)
}

# EDIT: "Control" here must match the prefix of your control_condition
annotation_col <- cd |>
  dplyr::mutate(
    condition = as.character(.data$condition),
    treatment = dplyr::if_else(stringr::str_detect(.data$condition, "^Control"), "Control",
                               stringr::str_extract(.data$condition, "^[A-Za-z]+")),
    time = dplyr::if_else(.data$treatment == "Control", "0",
                          stringr::str_extract(.data$condition, "(?<=_)[0-9]+"))
  ) |>
  dplyr::select(sample, treatment, time)

rownames(annotation_col) <- annotation_col$sample
annotation_col$sample <- NULL


# Order columns: Control first, then by time
# EDIT: replace with your own treatment names/order
treat_levels <- intersect(c("Control","TreatmentA","TreatmentC"), unique(annotation_col$treatment))

desired_order <- annotation_col |>
  dplyr::mutate(
    treatment = factor(.data$treatment, levels = treat_levels),
    time = as.numeric(.data$time)
  ) |>
  dplyr::arrange(.data$treatment, .data$time) |>
  rownames()

mat_scaled <- mat_scaled[, desired_order, drop = FALSE]
annotation_col <- annotation_col[desired_order, , drop = FALSE]

pheatmap::pheatmap(
  mat_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  show_colnames = FALSE
)



# ---- log2FC heatmap vs control (non-sig black), biomaRt gene labels ----
# Human dataset below. For mouse, use "mmusculus_gene_ensembl"
mart <- biomaRt::useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# EDIT: match the significance cutoffs used elsewhere
lfc_cutoff  <- 1.5
padj_cutoff <- 0.05
top_n_genes <- 250

# Contrasts vs control (set in script 1)
control_condition <- "Control_0"
keep_names <- names(res_shrunk_list)[stringr::str_detect(names(res_shrunk_list), paste0("_vs_", control_condition, "$"))]
res_use <- res_shrunk_list[keep_names]

# Convert results to data frames, keep contrast name
res_use_df <- purrr::imap(res_use, function(res, contrast_name) {
  df <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  df$contrast   <- contrast_name
  df
})

# Combine contrasts, keep genes meeting thresholds
sig_tbl <- dplyr::bind_rows(res_use_df) |>
  dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) |>
  dplyr::filter(padj < padj_cutoff, abs(log2FoldChange) >= lfc_cutoff) |>
  dplyr::select(ensembl_id, log2FoldChange, padj, contrast)

# Rank genes by max abs log2FC across contrasts
top_genes <- sig_tbl |>
  dplyr::group_by(ensembl_id) |>
  dplyr::summarise(max_abs_lfc = max(abs(log2FoldChange)), .groups = "drop") |>
  dplyr::arrange(dplyr::desc(max_abs_lfc)) |>
  dplyr::slice_head(n = top_n_genes) |>
  dplyr::pull(ensembl_id)

# Build log2FC matrix; non-significant values become NA
lfc_mat <- sapply(res_use_df, function(df) {
  idx <- match(top_genes, df$ensembl_id)
  lfc <- df$log2FoldChange[idx]
  p   <- df$padj[idx]
  
  is_sig <- !is.na(p) & !is.na(lfc) & (p < padj_cutoff) & (abs(lfc) >= lfc_cutoff)
  lfc[!is_sig] <- NA_real_
  
  lfc
})
colnames(lfc_mat) <- keep_names

# Parse contrast names into simplified labels (e.g. "TreatmentA 1")
contrast_tbl <- tibble::tibble(contrast = colnames(lfc_mat)) |>
  dplyr::mutate(
    condition = stringr::str_match(contrast, paste0("condition_([^_]+_[0-9]+)_vs_", control_condition))[, 2],
    treatment = stringr::str_extract(condition, "^[A-Za-z]+"),
    time      = as.numeric(stringr::str_extract(condition, "(?<=_)[0-9]+")),
    label     = paste(treatment, time)
  )

# Order columns by time then treatment. EDIT to match your treatments
treat_levels <- c("TreatmentA", "TreatmentC")
contrast_order <- contrast_tbl |>
  dplyr::mutate(
    time = factor(time, levels = sort(unique(time))),
    treatment = factor(treatment, levels = c(treat_levels, setdiff(unique(treatment), treat_levels)))
  ) |>
  dplyr::arrange(time, treatment) |>
  dplyr::pull(contrast)

lfc_mat <- lfc_mat[, contrast_order, drop = FALSE]

# Apply simplified column labels
col_label_map <- contrast_tbl$label
names(col_label_map) <- contrast_tbl$contrast
colnames(lfc_mat) <- col_label_map[colnames(lfc_mat)]

# Map Ensembl IDs to HGNC symbols, fall back to Ensembl ID
top_genes_clean <- sub("\\..*$", "", top_genes)

bm <- biomaRt::getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", "external_gene_name"),
  filters    = "ensembl_gene_id",
  values     = unique(top_genes_clean),
  mart       = mart
) |>
  tibble::as_tibble() |>
  dplyr::mutate(gene_name = dplyr::coalesce(hgnc_symbol, external_gene_name)) |>
  dplyr::filter(!is.na(gene_name) & gene_name != "") |>
  dplyr::distinct(ensembl_gene_id, .keep_all = TRUE)

gene_labels <- bm$gene_name[match(top_genes_clean, bm$ensembl_gene_id)]
gene_labels[is.na(gene_labels) | gene_labels == ""] <- top_genes[is.na(gene_labels) | gene_labels == ""]
rownames(lfc_mat) <- gene_labels

# Replace NA with zero only for clustering distance calculation
lfc_for_cluster <- lfc_mat
lfc_for_cluster[is.na(lfc_for_cluster)] <- 0
row_hc <- stats::hclust(stats::dist(lfc_for_cluster))

# Generate heatmap with non-significant values displayed in black
pheatmap::pheatmap(
  lfc_mat,
  cluster_rows = row_hc,
  cluster_cols = FALSE,
  na_col = "black",
  show_colnames = TRUE,
  show_rownames = FALSE,
  angle_col = 0,
  main = paste0("Top ", top_n_genes, " DE Genes Relative to Control (log2 Fold Change)")
)

# ---- Heatmap specific pathways ----

res_shrunk_list <- readRDS("projects/PROJECT_NAME/results/resShrink_annotated.rds")
# EDIT: match the significance cutoffs used elsewhere
lfc_cutoff  <- 1.5
padj_cutoff <- 0.05

# EDIT: pick a hallmark set from script 3's GSEA results
hallmark_name <- "<HALLMARK_SET_NAME>"
symbols <- unique(hallmark_sets[[hallmark_name]])

# Keep contrasts vs control
keep_names <- names(res_shrunk_list)[grepl(paste0("_vs_", control_condition, "$"), names(res_shrunk_list))]
res_use <- res_shrunk_list[keep_names]

# Strip Ensembl version suffix (ENSG...12 -> ENSG...)
clean_ens <- function(x) sub("\\..*$", "", as.character(x))

# Detect symbol/Ensembl ID columns
df0 <- as.data.frame(res_use[[1]])
sym_col <- intersect(c("SYMBOL", "symbol", "hgnc_symbol"), names(df0))[1]
ens_col <- intersect(c("ensembl_id", "ensembl_gene_id"), names(df0))[1]

# Convert to data frames, keep contrast name
res_use_df <- lapply(names(res_use), function(nm) {
  df <- as.data.frame(res_use[[nm]])
  df$contrast <- nm
  df
})
names(res_use_df) <- names(res_use)

# Genes in the hallmark set, significant in any contrast
sig_ens <- unique(unlist(lapply(res_use_df, function(df) {
  df$symbol <- df[[sym_col]]
  df$ensembl_id_clean <- clean_ens(df[[ens_col]])
  df <- df[!is.na(df$symbol) & df$symbol %in% symbols, ]
  df <- df[!is.na(df$padj) & !is.na(df$log2FoldChange), ]
  df <- df[df$padj < padj_cutoff & abs(df$log2FoldChange) >= lfc_cutoff, ]
  unique(df$ensembl_id_clean)
}), use.names = FALSE))

# Ensembl -> symbol map from the first contrast
sym_map <- res_use_df[[1]]
sym_map <- data.frame(
  ensembl_id_clean = clean_ens(sym_map[[ens_col]]),
  symbol = as.character(sym_map[[sym_col]]),
  stringsAsFactors = FALSE
)
sym_map <- sym_map[!is.na(sym_map$ensembl_id_clean) & !duplicated(sym_map$ensembl_id_clean), ]

# Normalized counts, with cleaned Ensembl IDs
norm_counts <- as.data.frame(DESeq2::counts(dds, normalized = TRUE))
norm_counts$ensembl_id_clean <- clean_ens(rownames(norm_counts))

# Keep only significant genes
norm_counts <- norm_counts[norm_counts$ensembl_id_clean %in% sig_ens, , drop = FALSE]
count_cols <- setdiff(names(norm_counts), "ensembl_id_clean")

# Wide to long: one row per gene-sample
counts_long <- data.frame(
  ensembl_id_clean = rep(norm_counts$ensembl_id_clean, times = length(count_cols)),
  sample = rep(count_cols, each = nrow(norm_counts)),
  norm_count = as.vector(as.matrix(norm_counts[, count_cols, drop = FALSE])),
  stringsAsFactors = FALSE
)

# Add sample metadata
coldata_df <- as.data.frame(SummarizedExperiment::colData(dds))
coldata_df$sample <- rownames(coldata_df)

# Combine into plotting table
plot_df <- dplyr::left_join(counts_long, coldata_df, by = "sample")
plot_df <- dplyr::left_join(plot_df, sym_map, by = "ensembl_id_clean")
plot_df <- plot_df[!is.na(plot_df$symbol) & plot_df$symbol %in% symbols, ]

# Parse condition into treatment and time (Control = time 0)
# EDIT: "Control" must match your control_condition prefix
plot_df <- plot_df |>
  dplyr::mutate(
    treatment = dplyr::if_else(grepl("^Control", condition), "Control", sub("_.*$", "", condition)),
    time = dplyr::if_else(grepl("^Control", condition), 0L, as.integer(sub("^.*_", "", condition))),
    group = paste0(treatment, "_", time)
  )

# Set x-axis ordering by time, then by treatment
group_order <- plot_df |>
  dplyr::distinct(group, treatment, time) |>
  dplyr::mutate(
    time = factor(time, levels = sort(unique(time))),
    treatment = factor(treatment, levels = c("Control", "TreatmentA", "TreatmentC"))
  ) |>
  dplyr::arrange(time, treatment) |>
  dplyr::pull(group)

plot_df <- plot_df |>
  dplyr::mutate(group = factor(group, levels = group_order))

# Mean/SEM per gene and group
summary_df <- plot_df |>
  dplyr::group_by(symbol, group, treatment) |>
  dplyr::summarise(
    n = sum(!is.na(norm_count)),
    mean_norm = mean(norm_count, na.rm = TRUE),
    sd_norm = stats::sd(norm_count, na.rm = TRUE),
    sem_norm = sd_norm / sqrt(n),
    .groups = "drop"
  )

# Color by treatment, or by group/time
color_mode <- "treatment"
file_ext   <- "png"

plot_df <- plot_df |>
  dplyr::mutate(color_key = if (color_mode == "group") as.character(group) else as.character(treatment))

summary_df <- summary_df |>
  dplyr::mutate(color_key = if (color_mode == "group") as.character(group) else as.character(treatment))

# Fill colors for bars/points
# EDIT: replace with your own treatment names/colors
fill_map <- if (color_mode == "treatment") {
  c(Control = "grey50", TreatmentA = "dodgerblue3", TreatmentC = "purple3")
} else {
  setNames(
    rep_len(c("grey50", "dodgerblue3", "purple3", "seagreen3", "orange3"), length(group_order)),
    group_order
  )
}

# Build a bar + jitter plot for a single gene symbol
make_gene_plot <- function(gene_symbol) {
  df_pts <- plot_df |> dplyr::filter(symbol == gene_symbol)
  df_sum <- summary_df |> dplyr::filter(symbol == gene_symbol)
  
  ggplot2::ggplot(df_sum, ggplot2::aes(x = group, y = mean_norm, fill = color_key)) +
    ggplot2::geom_col(width = 0.75) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_norm - sem_norm, ymax = mean_norm + sem_norm),
      width = 0.25
    ) +
    ggplot2::geom_point(
      data = df_pts,
      mapping = ggplot2::aes(x = group, y = norm_count, fill = color_key),
      position = ggplot2::position_jitter(width = 0.15, height = 0),
      shape = 21,
      size = 2.5,
      color = "black",
      stroke = 0.6
    ) +
    ggplot2::scale_fill_manual(values = fill_map) +
    ggplot2::labs(x = NULL, y = "Normalized counts", title = gene_symbol) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}

# Genes to plot, one figure each
sig_symbols <- summary_df |>
  dplyr::filter(!is.na(symbol) & symbol != "") |>
  dplyr::distinct(symbol) |>
  dplyr::pull(symbol)

# Output folder per hallmark set
outdir <- file.path("projects/PROJECT_NAME/figures", "bargraphs", hallmark_name)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Save one plot per gene symbol
for (g in sig_symbols) {
  p <- make_gene_plot(g)
  fn <- file.path(outdir, paste0(g, ".", file_ext))
  ggplot2::ggsave(fn, p, width = 8, height = 4, dpi = 300)
}




# ---- Venn Diagram ----
res_shrunk_list <- readRDS("projects/PROJECT_NAME/results/resShrink_annotated.rds")
# EDIT: match the significance cutoffs used elsewhere
lfc_cutoff  <- 1.5
padj_cutoff <- 0.05

# Contrasts vs control
keep_names <- names(res_shrunk_list)
keep_names <- keep_names[grepl(paste0("_vs_", control_condition, "$"), keep_names)]

# Parse treatment and time from contrast names
contrast_tbl <- data.frame(
  contrast  = keep_names,
  condition = sub(paste0("condition_([^_]+_[0-9]+)_vs_", control_condition), "\\1", keep_names),
  stringsAsFactors = FALSE
)

contrast_tbl$treatment <- sub("_.*", "", contrast_tbl$condition)
contrast_tbl$time      <- as.numeric(sub(".*_", "", contrast_tbl$condition))
# EDIT: compares TreatmentA vs TreatmentC. Adjust names/count for your arms.
# TreatmentB excluded since it was dropped during QC in script 1.
contrast_tbl <- contrast_tbl[contrast_tbl$treatment %in% c("TreatmentA", "TreatmentC"), ]

# Build significant gene sets per contrast
sig_sets <- list()

for (nm in keep_names) {
  res <- res_shrunk_list[[nm]]
  df  <- as.data.frame(res)
  df$ensembl_id <- rownames(df)
  
  sig_genes <- df$ensembl_id[
    !is.na(df$padj) &
      !is.na(df$log2FoldChange) &
      df$padj < padj_cutoff &
      abs(df$log2FoldChange) >= lfc_cutoff
  ]
  
  sig_sets[[nm]] <- unique(sub("\\..*$", "", sig_genes))
}

# Create Venn plots per timepoint
# EDIT: replace c(1, 2, 3) with your actual timepoint values
times <- sort(unique(contrast_tbl$time))
times <- times[times %in% c(1, 2, 3)]

venn_plots <- list()

for (tp in times) {

  a_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TreatmentA"]
  c_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TreatmentC"]

  sets_tp <- list(
    TreatmentA = sig_sets[[a_con]],
    TreatmentC = sig_sets[[c_con]]
  )

  venn_plots[[as.character(tp)]] <-
    ggvenn(
      sets_tp,
      show_elements = FALSE
    ) +
    ggtitle(paste0("T", tp)) +
    theme(plot.title = element_text(hjust = 0.5),
          size = 40)
}

# Draw all plots in a single row
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, length(venn_plots))))

for (i in seq_along(venn_plots)) {
  print(
    venn_plots[[i]],
    vp = viewport(layout.pos.row = 1, layout.pos.col = i)
  )
}


out_dir <- "projects/PROJECT_NAME/results"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

wb <- createWorkbook()

for (tp in times) {

  a_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TreatmentA"]
  c_con <- contrast_tbl$contrast[contrast_tbl$time == tp & contrast_tbl$treatment == "TreatmentC"]

  a_genes <- sig_sets[[a_con]]
  c_genes <- sig_sets[[c_con]]

  all_genes <- unique(c(a_genes, c_genes))

  a_df <- as.data.frame(res_shrunk_list[[a_con]])
  a_df$entrez_id <- sub("\\..*$", "", rownames(a_df))
  a_df <- a_df[!duplicated(a_df$entrez_id),
               c("entrez_id", "log2FoldChange", "padj")]
  names(a_df)[2:3] <- c("TreatmentA_log2FC", "TreatmentA_padj")

  c_df <- as.data.frame(res_shrunk_list[[c_con]])
  c_df$entrez_id <- sub("\\..*$", "", rownames(c_df))
  c_df <- c_df[!duplicated(c_df$entrez_id),
               c("entrez_id", "log2FoldChange", "padj")]
  names(c_df)[2:3] <- c("TreatmentC_log2FC", "TreatmentC_padj")

  out_df <- data.frame(entrez_id = all_genes, stringsAsFactors = FALSE)

  out_df$symbol_id <- mapIds(
    org.Hs.eg.db,
    keys      = as.character(out_df$entrez_id),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )

  out_df$membership <- ifelse(
    out_df$entrez_id %in% a_genes & out_df$entrez_id %in% c_genes, "Both",
    ifelse(out_df$entrez_id %in% a_genes, "TreatmentA_only", "TreatmentC_only")
  )

  out_df <- merge(out_df, a_df, by = "entrez_id", all.x = TRUE)
  out_df <- merge(out_df, c_df, by = "entrez_id", all.x = TRUE)

  out_df <- out_df[, c("entrez_id", "symbol_id", "membership",
                       "TreatmentA_log2FC", "TreatmentA_padj",
                       "TreatmentC_log2FC", "TreatmentC_padj")]

  out_df <- out_df[order(factor(out_df$membership,
                                levels = c("Both", "TreatmentA_only", "TreatmentC_only")),
                         out_df$entrez_id), ]

  sheet_name <- paste0("T", tp)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, out_df)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = seq_len(ncol(out_df)), widths = "auto")
}

saveWorkbook(wb, file.path(out_dir, "DE_genes_venn.xlsx"), overwrite = TRUE)