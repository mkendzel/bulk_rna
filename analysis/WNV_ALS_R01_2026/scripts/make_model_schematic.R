# Schematic of the limma-voom model, sized from the checkpoint dimensions so it
# cannot drift from the pipeline. Documentation figure, not part of the analysis
# flow - re-run after a design or filtering change.
# Run from the repo root.

# ---- Library ----
library(ggplot2)
library(dplyr)
library(tibble)

invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

EXP <- "spinal"   # the schematic is drawn for the spinal fit

# ---- Real numbers ----
counts <- load_checkpoint("annotated_counts", dir = dir_rds)
lc     <- load_checkpoint(paste0("logcpm_", EXP), dir = dir_rds)
sm     <- load_checkpoint(paste0("sample_meta_", EXP), dir = dir_rds)
vm     <- load_checkpoint(paste0("voom_", EXP), dir = dir_rds)
tt     <- load_checkpoint(paste0("tt_", EXP, "_annotated"), dir = dir_rds)

N_RAW      <- nrow(counts)
N_KEPT     <- nrow(lc)
N_SAMPLES  <- ncol(lc)
N_GROUPS   <- ncol(vm$design)
N_CONTRAST <- length(tt)

fmt <- function(x) format(x, big.mark = ",")

message(sprintf("%s: %s genes -> %s kept, %d samples, %d groups, %d contrasts",
                EXP, fmt(N_RAW), fmt(N_KEPT), N_SAMPLES, N_GROUPS, N_CONTRAST))

# ---- Drawing helpers ----
# Canvas is 0-100 in both directions; every position below is on that grid.
INK  <- "grey20"
MUTE <- "grey45"
LINE <- "grey55"

box <- function(x, y, w, h, fill = "white", colour = INK, lwd = 0.5, r = NULL) {
  annotate("rect", xmin = x, xmax = x + w, ymin = y, ymax = y + h,
           fill = fill, colour = colour, linewidth = lwd)
}

lab <- function(x, y, text, size = 9, face = "plain", colour = INK,
                hjust = 0.5, vjust = 0.5, family = "") {
  annotate("text", x = x, y = y, label = text, size = size / .pt,
           fontface = face, colour = colour, hjust = hjust, vjust = vjust,
           family = family, lineheight = 1.15)
}

arr <- function(x1, y1, x2, y2, colour = LINE, lwd = 0.6) {
  annotate("segment", x = x1, y = y1, xend = x2, yend = y2,
           arrow = arrow(length = unit(0.16, "cm"), type = "closed"),
           colour = colour, linewidth = lwd)
}

# ---- Row A: the pipeline ----
steps <- tibble::tribble(
  ~title,            ~detail,
  "counts",          sprintf("%s genes\n%d samples", fmt(N_RAW), N_SAMPLES),
  "filterByExpr\n+ TMM", sprintf("%s genes kept\nlibrary sizes scaled", fmt(N_KEPT)),
  "voom",            "log2 CPM\n+ precision weight\nper observation",
  "lmFit",           sprintf("~ 0 + condition\n%d group means", N_GROUPS),
  "contrasts.fit",   sprintf("%d subtractions\nbetween those means", N_CONTRAST),
  "eBayes",          "variance shrunk\nacross all genes",
  "topTable",        "logFC, t,\nP.Value, adj.P.Val"
)

NB <- nrow(steps)
BW <- 11.6; BH <- 13; GAP <- (100 - 8 - NB * BW) / (NB - 1)
bx <- 4 + (seq_len(NB) - 1) * (BW + GAP)
by <- 78

p <- ggplot() + xlim(0, 100) + ylim(0, 100)

# stage boxes; the three model stages are tinted with the contrast-type palette
stage_fill <- c("grey96", "grey96", "grey96",
                scales::alpha(type_cols[["within_line"]], 0.13),
                scales::alpha(type_cols[["genotype"]], 0.13),
                "grey96", "grey96")

for (i in seq_len(NB)) {
  p <- p +
    box(bx[i], by, BW, BH, fill = stage_fill[i]) +
    lab(bx[i] + BW / 2, by + BH - 3.1, steps$title[i], size = 10, face = "bold") +
    lab(bx[i] + BW / 2, by + BH / 2 - 2.2, steps$detail[i], size = 7.4, colour = MUTE)
  if (i < NB) p <- p + arr(bx[i] + BW + 0.6, by + BH / 2, bx[i + 1] - 0.6, by + BH / 2)
}

p <- p +
  lab(4, 96.5, "The limma-voom model", size = 15, face = "bold", hjust = 0) +
  lab(4, 93.2,
      sprintf("How %d spinal samples become %d contrasts — %s",
              N_SAMPLES, N_CONTRAST, "drawn from the actual checkpoint dimensions"),
      size = 9, colour = MUTE, hjust = 0)

# ---- Row B left: the cell-means design ----
# Shared vertical rhythm for both Row B panels: heading, then explanatory text
# hanging below it, then content starting at CONTENT_TOP.
HEAD_Y      <- 69
SUB_Y       <- 65.5
CONTENT_TOP <- 59

gx <- 4; gy <- 30; CW <- 9.4; CH <- 8.0
lines_use <- c("C9", "Ctrl3", "Ctrl2")
stims_use <- STIM_LEVELS

p <- p + lab(gx, HEAD_Y,
             sprintf("1.  Estimate one mean per group  (%d groups)", N_GROUPS),
             size = 11, face = "bold", hjust = 0) +
  lab(gx, SUB_Y,
      "The design has no intercept, so each cell below is a coefficient the model\nfits directly. Numbers are surviving replicates.",
      size = 8, colour = MUTE, hjust = 0, vjust = 1)

for (j in seq_along(stims_use))
  p <- p + lab(gx + (j - 0.5) * CW, gy + 3 * CH + 1.6, stims_use[j], size = 8.6, face = "bold")

for (i in seq_along(lines_use)) {
  yy <- gy + (3 - i) * CH
  p <- p + lab(gx - 1.2, yy + CH / 2, lines_use[i], size = 8.6, face = "bold",
               hjust = 1, colour = line_cols[[lines_use[i]]])
  for (j in seq_along(stims_use)) {
    cond <- paste(lines_use[i], stims_use[j], sep = "_")
    n    <- sum(sm$condition == cond)
    p <- p +
      box(gx + (j - 1) * CW, yy, CW - 0.6, CH - 0.6,
          fill = scales::alpha(line_cols[[lines_use[i]]], if (n == 1) 0.10 else 0.22),
          colour = line_cols[[lines_use[i]]], lwd = 0.45) +
      lab(gx + (j - 0.5) * CW - 0.3, yy + CH / 2 - 0.3,
          if (n == 0) "—" else paste0("n = ", n),
          size = 8.2, colour = if (n == 1) "#B22222" else INK,
          face = if (n == 1) "bold" else "plain")
  }
}

# ---- Row B right: contrasts as subtractions ----
cx <- 52
p <- p + lab(cx, HEAD_Y, "2.  Subtract those means", size = 11,
             face = "bold", hjust = 0) +
  lab(cx, SUB_Y,
      "Every contrast is a linear combination of the cells on the left.\nNothing is refitted; only the weights change.",
      size = 8, colour = MUTE, hjust = 0, vjust = 1)

fams <- tibble::tribble(
  ~type,          ~n,  ~q,                                  ~formula,
  "within_line",  9L,  "Did the stimulus do anything?",     "C9_IFNb  −  C9_Mock",
  "genotype",     8L,  "Do the lines differ, same stimulus?", "C9_IFNb  −  Ctrl3_IFNb",
  "interaction",  6L,  "Do the lines respond differently?",  "( C9_IFNb − C9_Mock )  −  ( Ctrl3_IFNb − Ctrl3_Mock )"
)

FH <- 10.5
for (i in seq_len(nrow(fams))) {
  yy  <- CONTENT_TOP - i * (FH + 1.2)
  col <- type_cols[[fams$type[i]]]
  p <- p +
    box(cx, yy, 45, FH, fill = scales::alpha(col, 0.10), colour = col, lwd = 0.45) +
    annotate("rect", xmin = cx, xmax = cx + 1.1, ymin = yy, ymax = yy + FH,
             fill = col, colour = NA) +
    lab(cx + 2.6, yy + FH - 2.9,
        sprintf("%s  (%d)", fams$type[i], fams$n[i]),
        size = 8.6, face = "bold", hjust = 0, colour = col) +
    lab(cx + 2.6, yy + FH - 5.9, fams$q[i], size = 7.6, hjust = 0, colour = MUTE) +
    lab(cx + 2.6, yy + 2.0, fams$formula[i], size = 7.4, hjust = 0,
        family = "mono", colour = INK)
}

# ---- finish ----

p <- p + theme_void() + theme(plot.margin = margin(2, 2, 2, 2))

save_fig(p, paste0("model_schematic_", EXP), width = 15, height = 9.5,
         subdir = "overview")
save_fig(p, paste0("model_schematic_", EXP), width = 15, height = 9.5,
         subdir = "overview", formats = FIG_FORMATS_PUB)

message("Wrote model_schematic_", EXP)
