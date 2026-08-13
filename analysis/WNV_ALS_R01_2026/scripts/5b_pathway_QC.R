# QC for the Figure 7 pathway panels in 5_Graphs.R. Reads the precomputed fgsea
# results, measures pathway-vs-pathway redundancy within each day and classifies
# pathway names by keyword. Flags only; nothing is dropped.
# Run from the repo root.

# ---- Library ----
library(fgsea)
library(dplyr)
library(tibble)
library(purrr)
library(stringr)
library(pheatmap)

# ---- Load helper functions ----
invisible(sapply(list.files("R", full.names = TRUE), source))
source("analysis/WNV_ALS_R01_2026/scripts/0_config.R")

# ---- Parameters ----
QC_OUT_DIR  <- file.path("outputs", "qc")
QC_GMT_PATH <- "genesets/h.all.v2025.1.Hs.symbols.gmt"

# Plotting cuts. Mirror FIG7_NES_CUT and FIG7_PVAL_CUT in 5_Graphs.R; keep in sync by hand.
QC_NES_CUT  <- 1.5
QC_PVAL_CUT <- 0.05

# hclust cut height on 1 - Jaccard of leading edge genes
QC_CUT_H <- 0.5

# A leading edge gene counts as promiscuous above this many other pathways
QC_PROMISC_N <- 1

# A pair is written to the overlap CSV when any metric clears its threshold
QC_PAIR_JAC  <- 0.10
QC_PAIR_CONT <- 0.20

# Leading edge genes listed per pathway, in fgsea's rank order
QC_TOP_LE <- 15

QC_SEED <- 1

# Pathway counts the current cuts produce. A mismatch means they drifted from 5_Graphs.R.
QC_EXPECTED_N <- c(d30 = 11L, d50 = 17L, d75 = 14L, d120 = 14L)

# No sampling or permutation; ordering is fixed by explicit tiebreaks.
set.seed(QC_SEED)

# ---- Keyword categories ----
# Scanned in order against the raw HALLMARK_* name, first match wins. Tissue terms
# precede generic ones so XENOBIOTIC_METABOLISM lands in hepatic_renal_gi, not metabolic.
QC_KEYWORDS <- tibble::tribble(
  ~pattern,                                                    ~category,                  ~flag,
  "MYOGENESIS|CARDIAC|MUSCLE",                                 "cardiac_muscle",           "review: tissue mismatch",
  "SPERMATOGENESIS|ANDROGEN|ESTROGEN",                         "reproductive",             "review: tissue mismatch",
  "PANCREAS|BILE_ACID|XENOBIOTIC",                             "hepatic_renal_gi",         "review: tissue mismatch",
  "COAGULATION|ANGIOGENESIS|HEME_METABOLISM",                  "hematopoietic_vascular",   "review: tissue mismatch",
  "EPITHELIAL_MESENCHYMAL|APICAL_JUNCTION|APICAL_SURFACE",     "epithelial_adhesion",      "review: tissue mismatch",
  "INTERFERON|INFLAMMATORY|TNFA|IL6|IL2|COMPLEMENT|ALLOGRAFT", "immune_inflammatory",      "review: cell-intrinsic vs infiltrate",
  "SYNAP|NEURO|AXON",                                          "neuronal_synaptic",        "",
  "GLIA|MYELIN|ASTRO|OLIGODEND",                               "glial",                    "",
  "HEDGEHOG|WNT_BETA_CATENIN|NOTCH|TGF_BETA",                  "developmental_patterning", "",
  "E2F|G2M|MITOTIC|MYC_TARGETS",                               "cell_cycle",               "",
  paste("P53|UNFOLDED_PROTEIN|REACTIVE_OXYGEN|UV_RESPONSE",
        "DNA_REPAIR|APOPTOSIS|HYPOXIA|PROTEIN_SECRETION", sep = "|"),
                                                               "stress_proteostasis",      "",
  paste("OXIDATIVE_PHOSPHORYLATION|GLYCOLYSIS|FATTY_ACID|CHOLESTEROL",
        "PEROXISOME|ADIPOGENESIS|MTORC1|PI3K", sep = "|"),
                                                               "metabolic",                ""
)

# First matching row of QC_KEYWORDS, with the substring that matched
classify_pathway <- function(pathway) {
  hit <- which(stringr::str_detect(pathway, QC_KEYWORDS$pattern))[1]

  if (is.na(hit)) {
    return(list(category = "other", matched_keyword = "", flag = "review: unclassified"))
  }

  list(
    category        = QC_KEYWORDS$category[hit],
    matched_keyword = stringr::str_extract(pathway, QC_KEYWORDS$pattern[hit]),
    flag            = QC_KEYWORDS$flag[hit]
  )
}

# ---- Data load ----
# Same directory, day map and file pattern as the Figure 7 block of 5_Graphs.R
fig7_dir  <- file.path(PROJ, "data", "fig7_als")
fig7_days <- c(d30 = "30", d50 = "50", d75 = "75", d120 = "120")

fig7_gsea <- setNames(lapply(names(fig7_days), function(d) {
  readRDS(file.path(fig7_dir, sprintf("fgsea_%s_ALS_vs_Ctrl(Batch).rds", d)))
}), names(fig7_days))

hallmark_sets <- gmtPathways(QC_GMT_PATH)

# ---- Plotted pathway set ----
# Same predicate as fig7_bubble_data(); its mutate() is cosmetic and adds no rows.
qc_plotted <- function(gsea_res) {
  gsea_res |>
    dplyr::filter(abs(NES) > QC_NES_CUT, pval < QC_PVAL_CUT) |>
    dplyr::arrange(NES) |>
    tibble::as_tibble()
}

# ---- Pairwise overlap ----
# Shared count, Jaccard and containment (|A n B| / min(|A|,|B|)) for every pair
overlap_metrics <- function(sets) {

  nm <- names(sets)
  n  <- length(nm)

  shared <- matrix(0L, n, n, dimnames = list(nm, nm))
  jac    <- matrix(0,  n, n, dimnames = list(nm, nm))
  cont   <- matrix(0,  n, n, dimnames = list(nm, nm))

  for (i in seq_len(n)) {
    for (j in seq_len(n)) {

      a <- sets[[i]]
      b <- sets[[j]]

      if (length(a) == 0 || length(b) == 0) next

      inter <- length(intersect(a, b))
      shared[i, j] <- inter
      jac[i, j]    <- inter / length(union(a, b))
      cont[i, j]   <- inter / min(length(a), length(b))
    }
  }

  list(shared = shared, jaccard = jac, containment = cont)
}

# Largest off-diagonal value per row, with the column it came from
max_partner <- function(m) {

  if (nrow(m) < 2) {
    return(list(value = setNames(rep(NA_real_, nrow(m)), rownames(m)),
                partner = setNames(rep(NA_character_, nrow(m)), rownames(m))))
  }

  off <- m
  diag(off) <- NA_real_

  list(
    value   = apply(off, 1, max, na.rm = TRUE),
    partner = colnames(off)[apply(off, 1, which.max)]
  )
}

# ---- Per-day QC ----
qc_day <- function(day) {

  k  <- qc_plotted(fig7_gsea[[day]])
  nm <- k$pathway
  n  <- length(nm)

  if (!is.na(QC_EXPECTED_N[day]) && n != QC_EXPECTED_N[day]) {
    warning(day, ": ", n, " pathways, expected ", QC_EXPECTED_N[day],
            " - cuts may have drifted from 5_Graphs.R", call. = FALSE)
  }

  le <- setNames(k$leadingEdge, nm)

  missing_sets <- setdiff(nm, names(hallmark_sets))
  if (length(missing_sets) > 0) {
    message(day, ": ", length(missing_sets), " pathway(s) absent from ", QC_GMT_PATH, ": ",
            paste(missing_sets, collapse = ", "))
  }
  full <- setNames(lapply(nm, function(p) if (is.null(hallmark_sets[[p]])) character(0)
                                          else hallmark_sets[[p]]), nm)

  ov_le  <- overlap_metrics(le)
  ov_set <- overlap_metrics(full)

  # Cluster on leading edge Jaccard
  if (n >= 2) {
    hc <- stats::hclust(stats::as.dist(1 - ov_le$jaccard), method = "average")
    cl <- stats::cutree(hc, h = QC_CUT_H)
  } else {
    hc <- NULL
    cl <- setNames(rep(1L, n), nm)
  }

  best_le  <- max_partner(ov_le$jaccard)
  best_set <- max_partner(ov_set$jaccard)

  # Plotted pathways containing each leading edge gene
  gene_n <- table(unlist(le, use.names = FALSE))

  promisc <- purrr::map_dfr(nm, function(p) {
    counts <- as.integer(gene_n[le[[p]]])
    tibble::tibble(
      pathway          = p,
      median_gene_n    = stats::median(counts),
      max_gene_n       = max(counts),
      frac_promiscuous = mean(counts - 1 > QC_PROMISC_N)
    )
  })

  cats <- purrr::map_dfr(nm, function(p) {
    cc <- classify_pathway(p)
    tibble::tibble(pathway = p, category = cc$category,
                   matched_keyword = cc$matched_keyword, relevance_flag = cc$flag)
  })

  qc_tbl <- k |>
    dplyr::transmute(
      pathway, size, NES, pval, padj,
      n_leading_edge = lengths(leadingEdge),
      cluster        = as.integer(cl[pathway])
    ) |>
    dplyr::group_by(cluster) |>
    dplyr::arrange(padj, dplyr::desc(abs(NES)), size, pathway, .by_group = TRUE) |>
    dplyr::mutate(
      is_representative = dplyr::row_number() == 1L,
      cluster_size      = dplyr::n()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      redundancy_flag = ifelse(cluster_size > 1L & !is_representative,
                               "redundancy candidate", ""),
      max_jaccard_le          = round(unname(best_le$value[pathway]), 3),
      max_jaccard_le_partner  = unname(best_le$partner[match(pathway, nm)]),
      max_jaccard_set         = round(unname(best_set$value[pathway]), 3),
      max_jaccard_set_partner = unname(best_set$partner[match(pathway, nm)])
    ) |>
    dplyr::left_join(promisc, by = "pathway") |>
    dplyr::left_join(cats, by = "pathway") |>
    dplyr::mutate(
      frac_promiscuous = round(frac_promiscuous, 3),
      top_leading_edge = vapply(pathway,
                                function(p) paste(utils::head(le[[p]], QC_TOP_LE), collapse = ";"),
                                character(1)),
      keep  = "",
      notes = ""
    ) |>
    dplyr::select(
      pathway, size, NES, pval, padj, n_leading_edge,
      cluster, cluster_size, is_representative, redundancy_flag,
      max_jaccard_le, max_jaccard_le_partner, max_jaccard_set, max_jaccard_set_partner,
      median_gene_n, max_gene_n, frac_promiscuous,
      category, matched_keyword, relevance_flag, top_leading_edge,
      keep, notes
    ) |>
    dplyr::arrange(cluster, dplyr::desc(is_representative), padj, pathway)

  # Long pair table, upper triangle only
  pairs <- if (n >= 2) {
    idx <- which(upper.tri(ov_le$jaccard), arr.ind = TRUE)

    tibble::tibble(
      day             = day,
      pathway_a       = nm[idx[, 1]],
      pathway_b       = nm[idx[, 2]],
      n_le_a          = lengths(le)[idx[, 1]],
      n_le_b          = lengths(le)[idx[, 2]],
      shared_le       = ov_le$shared[idx],
      jaccard_le      = round(ov_le$jaccard[idx], 3),
      containment_le  = round(ov_le$containment[idx], 3),
      size_a          = k$size[idx[, 1]],
      size_b          = k$size[idx[, 2]],
      shared_set      = ov_set$shared[idx],
      jaccard_set     = round(ov_set$jaccard[idx], 3),
      containment_set = round(ov_set$containment[idx], 3),
      same_cluster    = cl[nm[idx[, 1]]] == cl[nm[idx[, 2]]]
    ) |>
      dplyr::filter(jaccard_le  >= QC_PAIR_JAC  | containment_le  >= QC_PAIR_CONT |
                    jaccard_set >= QC_PAIR_JAC  | containment_set >= QC_PAIR_CONT) |>
      dplyr::arrange(dplyr::desc(jaccard_le), pathway_a, pathway_b)
  } else {
    tibble::tibble()
  }

  list(day = day, table = qc_tbl, pairs = pairs, jaccard_le = ov_le$jaccard,
       hc = hc, clusters = cl, n = n)
}

# ---- Run ----
ensure_dir(QC_OUT_DIR)

qc_all <- setNames(lapply(names(fig7_days), qc_day), names(fig7_days))

for (day in names(qc_all)) {

  res <- qc_all[[day]]

  utils::write.csv(res$table,
                   file.path(QC_OUT_DIR, sprintf("pathway_qc_table_%s.csv", day)),
                   row.names = FALSE)

  utils::write.csv(res$pairs,
                   file.path(QC_OUT_DIR, sprintf("pathway_pair_overlap_%s.csv", day)),
                   row.names = FALSE)

  if (res$n >= 2) {
    n_cl <- length(unique(res$clusters))
    pheatmap::pheatmap(
      res$jaccard_le,
      cluster_rows = res$hc,
      cluster_cols = res$hc,
      cutree_rows  = n_cl,
      cutree_cols  = n_cl,
      labels_row   = abbrev_pathway(rownames(res$jaccard_le)),
      labels_col   = abbrev_pathway(colnames(res$jaccard_le)),
      display_numbers = TRUE,
      number_format   = "%.2f",
      fontsize_number = 6,
      main     = sprintf("Leading edge Jaccard - Day %s (cut h = %s)", fig7_days[[day]], QC_CUT_H),
      filename = file.path(QC_OUT_DIR, sprintf("pathway_overlap_heatmap_%s.pdf", day)),
      width    = 9,
      height   = 8
    )
  }
}

# ---- Summary ----
for (day in names(qc_all)) {

  res <- qc_all[[day]]
  tb  <- res$table

  off <- res$jaccard_le
  if (res$n >= 2) diag(off) <- NA_real_

  message(sprintf(
    "%-5s %2d pathways | %2d clusters | %2d redundancy candidates | %2d tissue mismatch | max LE Jaccard %.2f",
    day, res$n, length(unique(res$clusters)),
    sum(tb$redundancy_flag == "redundancy candidate"),
    sum(tb$relevance_flag == "review: tissue mismatch"),
    if (res$n >= 2) max(off, na.rm = TRUE) else NA_real_
  ))
}

message("Wrote QC tables to ", QC_OUT_DIR)
