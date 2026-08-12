# Render the contrast registry to a one-page PDF reference sheet.
# Output: notes/WNV_ALS_R01_2026_contrasts.pdf
# Standalone utility - not part of the numbered pipeline. Re-run after editing
# contrast_registry in 0_config.R.

library(grid)
library(gridExtra)
library(dplyr)

source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

out_pdf <- file.path(PROJ, "notes", "WNV_ALS_R01_2026_contrasts.pdf")
ensure_dir(PROJ, "notes")

PAGE_W <- 11

reg <- contrast_registry |>
  dplyr::arrange(experiment, type, line, ref_line, stim)

# ---- Table styling ----
base_theme <- gridExtra::ttheme_default(
  core = list(
    fg_params = list(hjust = 0, x = 0.03, fontsize = 8, fontfamily = "mono"),
    padding   = unit(c(3, 2.4), "mm")
  ),
  colhead = list(
    fg_params = list(hjust = 0, x = 0.03, fontsize = 8.5, fontface = "bold"),
    bg_params = list(fill = "grey86"),
    padding   = unit(c(3, 2.4), "mm")
  )
)

# Shade unreplicated rows
themed_grob <- function(tbl, highlight) {
  th <- base_theme
  th$core$bg_params <- list(fill = ifelse(highlight, "#FBE4E4", "white"),
                            col = "grey80")
  gridExtra::tableGrob(tbl, rows = NULL, theme = th)
}

# Heading, subheading and a left-aligned table. Height is the exact sum of its
# parts so the table does not float in leftover space.
make_section <- function(heading, subheading, tbl) {

  tg <- themed_grob(tbl, tbl$n == 1)

  # pin the table to the left margin instead of centring it
  row <- gridExtra::arrangeGrob(
    tg, grid::nullGrob(),
    ncol   = 2,
    widths = unit.c(sum(tg$widths), unit(1, "null"))
  )

  heights <- unit.c(unit(2, "lines"), unit(1.15, "lines"), sum(tg$heights))

  grob <- gridExtra::arrangeGrob(
    grid::textGrob(heading, x = 0, hjust = 0,
                   gp = gpar(fontsize = 11.5, fontface = "bold")),
    grid::textGrob(subheading, x = 0, hjust = 0,
                   gp = gpar(fontsize = 8.5, col = "grey35", fontface = "italic")),
    row,
    ncol = 1, heights = heights
  )

  list(grob = grob, height = sum(heights))
}

fmt <- function(type_name) {
  reg |>
    dplyr::filter(type == type_name) |>
    dplyr::transmute(
      Contrast   = name,
      Arithmetic = gsub(" - ", " − ", expr),
      Exp        = experiment,
      n          = min_n
    ) |>
    as.data.frame()
}

w <- fmt("within_line")
g <- fmt("genotype")
i <- fmt("interaction")

n_by_line <- sample_map |> dplyr::count(line, name = "n")
n_of <- function(l) n_by_line$n[n_by_line$line == l]

subtitle <- paste0(
  "40 samples  |  spinal: C9 ", n_of("C9"), ", Ctrl3 ", n_of("Ctrl3"),
  ", Ctrl2 ", n_of("Ctrl2"), "    cortical: Cort ", n_of("Cort"),
  "  |  ", nrow(reg), " contrasts: ", nrow(w), " within-line, ",
  nrow(g), " genotype, ", nrow(i), " interaction",
  "  |  cutoffs adj.P.Val < ", P_CUT, ", |logFC| > ", LFC_CUT
)

footnote <- paste(
  "Batches are fit independently (separate DGEList, normalisation and dispersion). Cortical has no genotype or",
  "interaction contrasts: single line, no ALS counterpart, no WNV arm.",
  "Shaded rows are unreplicated (n = 1, from Ctrl2 IFNb/IFNg) - estimable via pooled variance and eBayes moderation,",
  "but with no within-group variance estimate. Treat as exploratory; the Ctrl3 rows are the replicated equivalents.",
  sep = "\n"
)

build <- function() {
  secs <- list(
    make_section("1.  Within-line - response to each stimulus",
                 "What does this stimulus do to this line?", w),
    make_section("2.  Genotype - ALS vs healthy at matched stimulus",
                 "How does the ALS line differ from a control line under the same condition?", g),
    make_section("3.  Interaction - differential response",
                 "Does the ALS line respond differently than the control? Baseline genotype differences cancel.", i)
  )

  heights <- unit.c(
    unit(2.2, "lines"), unit(1.6, "lines"),
    secs[[1]]$height, secs[[2]]$height, secs[[3]]$height,
    unit(5, "lines")
  )

  grobs <- list(
    grid::textGrob("WNV_ALS_R01_2026 - differential expression contrasts",
                   x = 0, hjust = 0, gp = gpar(fontsize = 15, fontface = "bold")),
    grid::textGrob(subtitle, x = 0, hjust = 0,
                   gp = gpar(fontsize = 8.5, col = "grey30")),
    secs[[1]]$grob, secs[[2]]$grob, secs[[3]]$grob,
    grid::textGrob(footnote, x = 0, hjust = 0,
                   gp = gpar(fontsize = 7.5, col = "grey35", lineheight = 1.35))
  )

  list(grobs = grobs, heights = heights)
}

# Measure on a throwaway device so the page is exactly as tall as its content
pdf(NULL, width = PAGE_W, height = 30)
layout <- build()
page_h <- grid::convertHeight(sum(layout$heights), "in", valueOnly = TRUE) + 0.7
dev.off()

pdf(out_pdf, width = PAGE_W, height = page_h)
layout <- build()
grid.arrange(grobs = layout$grobs, ncol = 1, heights = layout$heights,
             padding = unit(0, "line"),
             vp = grid::viewport(width = unit(1, "npc") - unit(0.7, "in"),
                                 height = unit(1, "npc") - unit(0.5, "in")))
dev.off()

cat("Wrote ", out_pdf, "  (", PAGE_W, " x ", round(page_h, 2), " in)\n", sep = "")
