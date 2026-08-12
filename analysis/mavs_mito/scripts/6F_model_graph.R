library(readxl)
library(dplyr)
library(tidyr)

path   <- "projects/mavs_mito/data/raw/AV_MAVS_SiT_p14.xlsx"
sheets <- c("Spleen", "Brain", "Spinal Cord")

# Read each tissue sheet, name the unnamed first column, tag the tissue, stack
raw <- lapply(sheets, function(s) {
  read_excel(path, sheet = s) |>
    rename(rep = 1) |>
    mutate(tissue = s)
}) |>
  bind_rows()

# Wide genotype columns to long: one row per mouse-tissue observation
long <- raw |>
  pivot_longer(c("WT", "MAVS-/-"), names_to = "genotype", values_to = "fc") |>
  filter(!is.na(fc)) |>
  mutate(
    genotype = factor(genotype, levels = c("WT", "MAVS-/-")),
    tissue   = factor(tissue, levels = sheets),
    log10_fc = log10(fc)
  )

# Welch two-sample t-test on log10 fold-change, WT vs MAVS-/-, within each tissue
welch <- lapply(sheets, function(tis) {
  d        <- subset(long, tissue == tis)
  tt       <- t.test(log10_fc ~ genotype, data = d)   # Welch, unequal variance
  diff_log <- -diff(tt$estimate)                      # mean(WT) - mean(MAVS-/-), log10 scale
  ci_log   <- tt$conf.int                             # CI for WT - MAVS-/-, log10 scale
  data.frame(
    tissue          = tis,
    log10_WT        = unname(tt$estimate[1]),
    log10_MAVS      = unname(tt$estimate[2]),
    fold_MAVS_vs_WT = 10^(-diff_log),                 # geometric-mean fold difference
    ci_low          = 10^(-ci_log[2]),
    ci_high         = 10^(-ci_log[1]),
    t               = unname(tt$statistic),
    df              = unname(tt$parameter),
    p.value         = tt$p.value
  )
}) |>
  bind_rows()

welch$p.holm <- p.adjust(welch$p.value, method = "holm")
welch