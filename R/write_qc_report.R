# Markdown QC report writer. Same accumulator pattern as render_checkpoint_md():
# build a character vector, writeLines it once.

# Format a metric value for display. compact = TRUE drops separators, suffixes and
# the percent sign so a value fits inside a plot tile.
fmt_qc <- function(value, unit, compact = FALSE) {
  num <- if (compact) {
    ifelse(abs(value) >= 1e6, sprintf("%.1fM", value / 1e6),
           ifelse(abs(value) >= 1e3, sprintf("%.1fK", value / 1e3),
                  formatC(round(value), format = "d")))
  } else {
    formatC(round(value), format = "d", big.mark = ",")
  }

  # rRNA runs four orders of magnitude below one percent, where %.2f prints 0.00
  pct <- if (compact) {
    sprintf("%.1f", value)
  } else {
    ifelse(value != 0 & abs(value) < 0.01,
           sprintf("%.2e%%", value),
           sprintf("%.2f%%", value))
  }

  ratio <- if (compact) sprintf("%.2f", value) else sprintf("%.3f", value)

  ifelse(is.na(value), "-",
         ifelse(unit == "pct",   pct,
                ifelse(unit == "ratio", ratio, num)))
}

# Report layout. One row per section, in output order. `category` selects rows of
# qc_status; NA means the section is plots only. `plots` is a comma-separated list
# of file names in the figure directory; NA means no plot.
QC_REPORT_SECTIONS <- tibble::tribble(
  ~category,     ~heading,                 ~plots,
  "libsize",     "Total counts",           "qc_total_counts.png",
  "complexity",  "Library complexity",     "qc_library_complexity.png",
  "detection",   "Gene detection",         "qc_gene_detection.png",
  "biotype",     "Gene biotypes",          "qc_biotypes.png",
  "alignment",   "Alignment",              NA_character_,
  NA_character_, "Variance stabilization", "qc_mean_sd.png",
  "clustering",  "Sample clustering",      "qc_sample_clustering.png",
  "pca",         "PCA",                    "qc_pca.png,qc_pca_scatters.png",
  "replicate",   "Replicate variability",  NA_character_
)

# Value with its status appended unless it passed, as used in every per-sample table
.qc_cell <- function(value, unit, status) {
  paste0(fmt_qc(value, unit),
         ifelse(status == "PASS" | status == "NA", "", paste0(" (", status, ")")))
}

# One sample-by-metric table for the given rows of qc_status
.qc_wide <- function(status_rows) {
  status_rows |>
    dplyr::mutate(cell = .qc_cell(value, unit, status)) |>
    dplyr::select(sample, metric, cell) |>
    tidyr::pivot_wider(names_from = metric, values_from = cell) |>
    dplyr::arrange(sample)
}

write_qc_report <- function(path,
                            qc_summary,
                            qc_status,
                            thresholds,
                            params,
                            prov,
                            fig_rel           = "../figures/qc",
                            counts_checkpoint = NA_character_,
                            sections          = QC_REPORT_SECTIONS) {

  L   <- character(0)
  add <- function(...) L[[length(L) + 1]] <<- paste0(...)
  tbl <- function(df) {
    for (line in knitr::kable(as.data.frame(df), format = "markdown")) add(line)
    add("")
  }

  n_status <- table(qc_status$status)

  add("# QC report — WNV_ALS_R01_2026")
  add("")
  tbl(data.frame(
    field = c("Generated", "Git", "R", "Platform", "Counts checkpoint",
              "Samples", "Checks", "Pass", "Warn", "Fail"),
    value = c(prov$timestamp,
              ifelse(is.na(prov$git), "-", prov$git),
              prov$r_version,
              prov$platform,
              ifelse(is.na(counts_checkpoint), "-", counts_checkpoint),
              dplyr::n_distinct(qc_status$sample),
              nrow(thresholds),
              n_status[["PASS"]], n_status[["WARN"]], n_status[["FAIL"]])
  ))

  add("## Status summary")
  add("")
  tbl(qc_status |>
        dplyr::group_by(metric, unit) |>
        dplyr::summarise(
          min    = fmt_qc(min(value, na.rm = TRUE), unit[1]),
          median = fmt_qc(stats::median(value, na.rm = TRUE), unit[1]),
          max    = fmt_qc(max(value, na.rm = TRUE), unit[1]),
          warn   = fmt_qc(warn_at[1], unit[1]),
          fail   = fmt_qc(fail_at[1], unit[1]),
          n_warn = sum(status == "WARN"),
          n_fail = sum(status == "FAIL"),
          .groups = "drop"
        ) |>
        dplyr::select(-unit))

  add("## Flagged samples")
  add("")
  flagged <- qc_status |>
    dplyr::filter(status %in% c("WARN", "FAIL")) |>
    dplyr::arrange(dplyr::desc(status), metric, value) |>
    dplyr::transmute(
      sample, metric, status,
      value = fmt_qc(value, unit),
      limit = paste0(ifelse(direction == "higher", "min ", "max "),
                     fmt_qc(ifelse(status == "FAIL", fail_at, warn_at), unit))
    )
  if (nrow(flagged) > 0) tbl(flagged) else { add("None."); add("") }

  add("## Status grid")
  add("")
  add("![](", fig_rel, "/qc_status_grid.png)")
  add("")

  add("## Per-sample metrics")
  add("")
  add("![](", fig_rel, "/qc_metrics_by_sample.png)")
  add("")

  for (i in seq_len(nrow(sections))) {
    sec <- sections[i, ]
    add("## ", sec$heading)
    add("")

    if (!is.na(sec$category)) {
      rows <- qc_status[qc_status$category == sec$category, ]
      rows <- droplevels(rows)
      if (nrow(rows) > 0) tbl(.qc_wide(rows))
    }

    if (!is.na(sec$plots)) {
      for (p in strsplit(sec$plots, ",")[[1]]) add("![](", fig_rel, "/", trimws(p), ")")
      add("")
    }
  }

  add("## Thresholds")
  add("")
  tbl(thresholds)

  add("## Parameters")
  add("")
  tbl(data.frame(parameter = names(params), value = unlist(params)))

  add("## Full metric table")
  add("")
  tbl(.qc_wide(qc_status))

  writeLines(unlist(L), path)
  invisible(path)
}
