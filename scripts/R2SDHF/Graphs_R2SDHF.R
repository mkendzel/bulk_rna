# ---- Triple Volcano: IFNa / STAT3 / OxPhos at TX48 (R2SDHF) ----

library(tidyverse)
library(fgsea)
library(ggrepel)
library(patchwork)

contrast <- "TX48"

base_dir     <- "data/Anderson-Suthar_collab/r_objects"
gsea_results <- readRDS(file.path(base_dir, paste0("fgsea_", contrast, "_vs_Mock(R2SDHF).rds")))
res_results  <- readRDS(file.path(base_dir, paste0("res_", contrast, "(R2SDHF).rds")))

time_display <- paste0(str_extract(contrast, "\\d+"), "h post Infection")

# Define panels with per-panel max_overlaps
panels <- list(
  list(regex = "INTERFERON_ALPHA",           title = "IFN-\u03b1 Response",        max_overlaps = 19),
  list(regex = "IL6_JAK_STAT3",              title = "IL-6/JAK/STAT3 Signaling",   max_overlaps = 21),
  list(regex = "OXIDATIVE_PHOSPHORYLATION",  title = "Oxidative Phosphorylation",  max_overlaps = 45)
)

# Base volcano data frame
volcano_df <- res_results %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    FC      = 2^logFC,
    log10_p = log10(P.Value)
  )

# Build one volcano plot per panel
make_volcano <- function(panel, volcano_df, gsea_results, time_display) {
  target_paths <- gsea_results %>%
    filter(str_detect(pathway, regex(panel$regex, ignore_case = TRUE)))
  target_genes <- unique(unlist(target_paths$leadingEdge))
  
  vdf <- volcano_df %>%
    mutate(
      is_target = gene %in% target_genes,
      target_direction = case_when(
        is_target & logFC > 0 ~ "Up",
        is_target & logFC < 0 ~ "Down",
        TRUE ~ NA_character_
      )
    )
  
  # Symmetric x limits based on max absolute fold change, capped at 32
  max_fc <- vdf %>%
    filter(is_target) %>%
    pull(FC) %>%
    {max(abs(log2(.)), na.rm = TRUE)} %>%
    min(5)
  x_lim <- c(2^(-max_fc) / 1.2, 2^(max_fc) * 1.2)
  
  # Dynamic x-axis breaks based on data range
  min_pow <- ceiling(log2(x_lim[1]))
  max_pow <- floor(log2(x_lim[2]))
  x_breaks <- 2^(min_pow:max_pow)
  x_labels <- sapply(x_breaks, function(x) {
    if (x == 1) "1"
    else if (x > 1) paste0("+", x)
    else paste0("-", round(1/x))
  })
  
  # Set y-axis floor based on the least significant target gene
  target_min_p <- vdf %>%
    filter(is_target) %>%
    pull(log10_p) %>%
    min(na.rm = TRUE)
  y_floor <- target_min_p * 1.05
  
  # Label all target genes, let ggrepel manage overlap
  label_df <- vdf %>%
    filter(is_target)
  
  ggplot(vdf, aes(x = FC, y = log10_p)) +
    geom_point(data = filter(vdf, !is_target),
               color = "grey70", size = 0.6, alpha = 0.4) +
    geom_point(data = filter(vdf, is_target),
               aes(color = target_direction), size = 2.2, alpha = 0.9) +
    geom_text_repel(data = label_df,
                    aes(label = gene), size = 3.5, fontface = "italic",
                    max.overlaps = panel$max_overlaps, segment.size = 0.3,
                    min.segment.length = 0, box.padding = 0.4,
                    force = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "black"), guide = "none") +
    scale_x_continuous(
      trans = "log2",
      breaks = x_breaks,
      labels = x_labels
    ) +
    coord_cartesian(xlim = x_lim, ylim = c(y_floor, 0)) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(x = "FC: WNV vs. Mock", y = "log10(p-value)",
         title = paste0(panel$title, " (", time_display, ")")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank()
    )
}

plot_list <- map(panels, make_volcano,
                 volcano_df = volcano_df,
                 gsea_results = gsea_results,
                 time_display = time_display)

# Combine three panels in a single row
combined <- wrap_plots(plot_list, nrow = 1) +
  plot_annotation(tag_levels = "A")

combined
ggsave("figures/Shilu_SERVC/WNV48_triple_volcano.png", combined,
       width = 16.66, height = 3.83, dpi = 300)

# ---- Triple Volcano: TX48 leading edge genes plotted on ROV48 (R2SDHF) ----

library(tidyverse)
library(fgsea)
library(ggrepel)
library(patchwork)

base_dir <- "data/Anderson-Suthar_collab/r_objects"

# Load TX48 GSEA to get leading edge genes
gsea_tx48 <- readRDS(file.path(base_dir, "fgsea_TX48_vs_Mock(R2SDHF).rds"))

# Load ROV48 DE results for the volcano background
res_rov48 <- readRDS(file.path(base_dir, "res_ROV48(R2SDHF).rds"))

time_display <- "48h post Infection"

# Define panels with per-panel max_overlaps
panels <- list(
  list(regex = "INTERFERON_ALPHA",           title = "IFN-\u03b1 Response",        max_overlaps = 19),
  list(regex = "IL6_JAK_STAT3",              title = "IL-6/JAK/STAT3 Signaling",   max_overlaps = 21),
  list(regex = "OXIDATIVE_PHOSPHORYLATION",  title = "Oxidative Phosphorylation",  max_overlaps = 45)
)

# Base volcano data frame from ROV48
volcano_df <- res_rov48 %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  mutate(
    FC      = 2^logFC,
    log10_p = log10(P.Value)
  )

# Build one volcano plot per panel
make_volcano <- function(panel, volcano_df, gsea_tx48, time_display) {
  # Pull leading edge genes from TX48 GSEA
  target_paths <- gsea_tx48 %>%
    filter(str_detect(pathway, regex(panel$regex, ignore_case = TRUE)))
  target_genes <- unique(unlist(target_paths$leadingEdge))
  
  # Flag those genes on the ROV48 volcano
  vdf <- volcano_df %>%
    mutate(
      is_target = gene %in% target_genes,
      target_direction = case_when(
        is_target & logFC > 0 ~ "Up",
        is_target & logFC < 0 ~ "Down",
        TRUE ~ NA_character_
      )
    )
  
  # Symmetric x limits based on max absolute fold change, capped at 32
  max_fc <- vdf %>%
    filter(is_target) %>%
    pull(FC) %>%
    {max(abs(log2(.)), na.rm = TRUE)} %>%
    min(5)
  x_lim <- c(2^(-max_fc) / 1.2, 2^(max_fc) * 1.2)
  
  # Dynamic x-axis breaks based on data range
  min_pow <- ceiling(log2(x_lim[1]))
  max_pow <- floor(log2(x_lim[2]))
  x_breaks <- 2^(min_pow:max_pow)
  x_labels <- sapply(x_breaks, function(x) {
    if (x == 1) "1"
    else if (x > 1) paste0("+", x)
    else paste0("-", round(1/x))
  })
  
  # Set y-axis floor based on the least significant target gene
  target_min_p <- vdf %>%
    filter(is_target) %>%
    pull(log10_p) %>%
    min(na.rm = TRUE)
  y_floor <- target_min_p * 1.05
  
  # Label all target genes, let ggrepel manage overlap
  label_df <- vdf %>%
    filter(is_target)
  
  ggplot(vdf, aes(x = FC, y = log10_p)) +
    geom_point(data = filter(vdf, !is_target),
               color = "grey70", size = 0.6, alpha = 0.4) +
    geom_point(data = filter(vdf, is_target),
               aes(color = target_direction), size = 2.2, alpha = 0.9) +
    geom_text_repel(data = label_df,
                    aes(label = gene), size = 3.5, fontface = "italic",
                    max.overlaps = panel$max_overlaps, segment.size = 0.3,
                    min.segment.length = 0, box.padding = 0.4,
                    force = 2) +
    scale_color_manual(values = c("Up" = "red", "Down" = "black"), guide = "none") +
    scale_x_continuous(
      trans = "log2",
      breaks = x_breaks,
      labels = x_labels
    ) +
    coord_cartesian(xlim = x_lim, ylim = c(y_floor, 0)) +
    geom_vline(xintercept = 1, color = "black", linewidth = 0.5) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(x = "FC: ROV vs. Mock", y = "log10(p-value)",
         title = paste0(panel$title, " (", time_display, ")")) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      panel.grid.minor = element_blank()
    )
}

plot_list <- map(panels, make_volcano,
                 volcano_df = volcano_df,
                 gsea_tx48 = gsea_tx48,
                 time_display = time_display)

# Combine three panels in a single row
combined <- wrap_plots(plot_list, nrow = 1) +
  plot_annotation(tag_levels = "A")

combined

ggsave("figures/Shilu_SERVC/ROV48_tx48genes_triple_volcano.png", combined,
       width = 16.66, height = 3.83, dpi = 300)