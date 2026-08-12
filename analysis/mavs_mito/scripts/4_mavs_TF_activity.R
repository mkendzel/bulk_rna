# =============================================
# 4 - Transcription factor activity (decoupleR)
# =============================================
# Infers TF activity from the brain WT vs MAVSKO D7 contrast using
# CollecTRI mouse regulons and a univariate linear model (run_ulm) on the
# moderated t statistic. Produces a signed ranking of the top TFs.
# Positive score = higher activity in WT; negative = higher in Mavs-/-.

# ---- Libraries ----
library(tidyverse)
library(decoupleR)
library(ggplot2)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))

# CollecTRI mouse regulons, complexes kept intact
net <- decoupleR::get_collectri(organism = "mouse", split_complexes = FALSE)

res_voom_brain_d0_d7_annotated <- load_checkpoint("res_voom_brain_d0_d7_annotated",
  dir = "projects/mavs_mito/data/r_objects")

contrast_tbl <- res_voom_brain_d0_d7_annotated$WT_vs_MAVSKO_D7

# Collapse to one row per symbol, keeping the most highly expressed duplicate
contrast_tbl <- contrast_tbl[!is.na(contrast_tbl$SYMBOL), ]
contrast_tbl <- contrast_tbl[order(-contrast_tbl$AveExpr), ]
contrast_tbl <- contrast_tbl[!duplicated(contrast_tbl$SYMBOL), ]

# Gene level statistic matrix using the moderated t as the input signal
mat <- matrix(contrast_tbl$t, ncol = 1,
              dimnames = list(contrast_tbl$SYMBOL, "WT_vs_MAVSKO"))

# Univariate linear model TF activity inference across the full transcriptome
tf_acts <- decoupleR::run_ulm(
  mat     = mat,
  net     = net,
  .source = "source",
  .target = "target",
  .mor    = "mor",
  minsize = 5
)

# Multiple testing correction and ranking
# Positive score is higher activity in WT, negative is higher in Mavs-/-
tf_acts <- as.data.frame(tf_acts)
tf_acts$padj <- p.adjust(tf_acts$p_value, method = "BH")
tf_acts <- tf_acts[order(-tf_acts$score), ]

# Top TFs by absolute activity for a signed ranking
n_show <- 20
top_tfs <- tf_acts[order(-abs(tf_acts$score)), ]
top_tfs <- top_tfs[seq_len(min(n_show, nrow(top_tfs))), ]
top_tfs$source <- factor(top_tfs$source, levels = top_tfs$source[order(top_tfs$score)])
top_tfs$direction <- ifelse(top_tfs$score > 0, "WT", "Mavs-/-")

p <- ggplot(top_tfs, aes(x = score, y = source, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c("WT" = "grey50", "Mavs-/-" = "royalblue3")) +
  labs(x = "TF activity (ULM score)", y = NULL, fill = NULL) +
  theme_bw()

p
