# WNV_ALS_R01_2026

Human brain organoid bulk RNA-seq: ALS (C9orf72) vs healthy control lines, challenged
with IFN-beta, IFN-gamma, or West Nile virus.

Run everything with the **repo root** (`bulk_rna/`) as the working directory — all paths
are relative to it, and scripts source shared helpers from `R/`.

## Design

40 samples, two batches modelled independently.

| Line | Sheet name | Experiment | Mock | IFNb | IFNg | WNV | n |
|---|---|---|---|---|---|---|---|
| `C9` | C9_1 ALS spinal organoids | spinal | 3 | 3 | 3 | 3 | 12 |
| `Ctrl3` | Ctrl_3 healthy control | spinal | 2 | 3 | 3 | 3 | 11 |
| `Ctrl2` | Ctrl_2 healthy control | spinal | 3 | **1** | **1** | 3 | 8 |
| `Cort` | Cortical organoids, uninfected | cort | 3 | 3 | 3 | — | 9 |

The cortical set is a separate experiment (`notes/hBO list for bulk RNA seq_SM_20260727.docx`:
"a separate set of experiment … PAIRWISE — all by itself"), so it gets its own `DGEList`,
normalisation factors and dispersion estimates. Nothing is pooled across batches.

**Ctrl2 IFNb and IFNg are n = 1.** The cell-means fit still estimates those group means
(residual variance is pooled across all 12 groups and moderated by `eBayes`), but the six
contrasts touching them are unreplicated. They carry `min_n = 1` in the contrast registry
and are drawn as triangles in the GSEA bubble plots. Treat them as exploratory.

## Sample naming

Samples were submitted to Plasmidsaurus under coded names. `sample_map` in
[`0_config.R`](0_config.R) is the authoritative translation, transcribed from
`notes/SM072026_bulk RNAseq.xlsx`.

```
plasmid grammar   <LINE>-<STIM><rep>_<globalIndex>      e.g. C9-TX1_10
line   C9 | CT3 -> Ctrl3 | CT2 -> Ctrl2 | CO -> Cort
stim   M -> Mock | b -> IFNb | g -> IFNg | TX -> WNV
result <line>_<stim>_<rep>                              e.g. C9_WNV_1
```

Line tokens contain no underscore, so `condition` is recoverable as
`sub("_[0-9]+$", "", sample)` and `_`-delimited splits stay unambiguous.

Plasmidsaurus does not ship the coded names in the matrix. Its columns are
`FCMSVB_<idx>_cpm` / `FCMSVB_<idx>_count`, where `<idx>` is the trailing global index of
the plasmid code — `FCMSVB_10` is `C9-TX1_10` is `C9_WNV_1`.

That identity is verified, not assumed. `vendor_key()` in [`0_config.R`](0_config.R) reads
the vendor's own idx → code map out of the per-sample filenames inside
`data/raw/FCMSVB_results.zip` (`FCMSVB_mapping-stats/FCMSVB_<idx>_<plasmid>.tsv`, listed
without extracting) and script 1 asserts it matches `sample_map`. If a future vendor run
numbers its columns independently of the notebook index, that assertion fails — without it
the import would complete with every sample label permuted and every other check still
passing.

Script 1 then matches columns in three tiers: vendor key, exact plasmid code, trailing
index. Unmapped columns are listed and then `stopifnot` fails rather than silently
dropping samples.

## Raw files

Only two files in `data/raw/` are used:

| File | Used for |
|---|---|
| `FCMSVB-expression-matrix.tsv` | counts and CPM, 40 samples |
| `FCMSVB_results.zip` | sample-name map, alignment read stats, MultiQC / PCA / correlation reports |

`FCMSVB_bam.zip` (14 GB) and `FCMSVB_fastq.zip` (31 GB) are **not used** — the vendor
already produced the count matrix and nothing here re-quantifies. Download them only if
raw reads are needed for a GEO/SRA deposition.

The results zip is read in place with `unzip(list = TRUE)` and `unz()`; nothing is
extracted into the repo. All of `data/raw/` is gitignored.

## QC

[`2_QC.R`](2_QC.R) runs 14 per-sample checks and writes `results/qc_report.md` plus ten
plots in `figures/qc/`. The report is generated on every run — tables and figures only, no
commentary. QC filters nothing and writes no checkpoint, so it can be re-run freely while
tuning thresholds; `drop_samples` lives in [`3_DE_limma.R`](3_DE_limma.R).

Read `qc_status_grid.png` for every sample against every check at a glance, and
`qc_metrics_by_sample.png` for the same data as bars with the cutoffs drawn in. The
RNAseqQC figures beside them describe the run as a whole rather than naming samples.

Checks come in three groups: vendor alignment stats read out of `FCMSVB_results.zip`, count
matrix statistics (library size, complexity, gene detection, biotype composition), and
vst-derived agreement measures (within-condition correlation, PCA centroid distance,
replicate deviation). Plot categories follow the
[RNAseqQC vignette](https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html);
`plot_chromosome()` is not used because the vendor matrix carries no gene coordinates.

**Every cutoff lives in the `QC_THRESHOLDS` tribble at the top of `2_QC.R`**, with
non-threshold knobs in `QC_PARAMS` beside it. All values are absolute. Nothing else in the
script hardcodes a limit.

The vst-derived cutoffs (`within_condition_cor`, `pca_centroid_dist`, `median_abs_M`) are
calibrated against this dataset's distribution and will need retuning on a different run.
`mapping_rate`, `dedup_rate`, `protein_coding_pct` and `rRNA_pct` are sentinels — they sit
far from anything these samples produce and are there to catch a future failure, not to
flag this one.

As of the 2026-08 run: 506 pass, 30 warn, 20 fail. Seven samples fail at least one check —
`Cort_Mock_1`, `Cort_Mock_2`, `Cort_Mock_3`, `Cort_IFNg_2`, `Ctrl3_IFNb_1`, `Ctrl3_IFNb_2`,
`Ctrl2_WNV_3`. `drop_samples` is empty.

## Contrasts

25 contrasts, built in `contrast_registry` ([`0_config.R`](0_config.R)). Scripts 4 and 5
join against that table on `name` to recover grouping metadata and plot labels.

| type | n | pattern |
|---|---|---|
| `within_line` | 11 | `<line>_<stim> - <line>_Mock`, all four lines |
| `genotype` | 8 | `C9_<stim> - Ctrl3_<stim>` and `C9_<stim> - Ctrl2_<stim>` |
| `interaction` | 6 | `(C9_<stim> - C9_Mock) - (<ref>_<stim> - <ref>_Mock)` |

The interaction family is why this pipeline uses **limma-voom** rather than DESeq2:
`makeContrasts` on a `~ 0 + condition` cell-means design expresses a difference of
differences directly, whereas `lfcShrink(type = "apeglm")` accepts only a single
coefficient.

Not built: cortical-vs-spinal mock comparison. The docx asks for it ("new batch to compare
with the first experiment mocks") but it confounds tissue with batch. Adding it is one more
registry row plus a combined fit.

## Scripts

| | |
|---|---|
| [`0_config.R`](0_config.R) | Paths, cutoffs, `sample_map`, `contrast_registry`, palettes. Sourced by 1–5. |
| [`1_df_count.R`](1_df_count.R) | Import → rename → checkpoint. Run once per delivery. |
| [`2_QC.R`](2_QC.R) | 14 QC checks → `results/qc_report.md` + `figures/qc/`. Writes no checkpoint. |
| [`3_DE_limma.R`](3_DE_limma.R) | `drop_samples` → split by experiment → voom/lmFit/eBayes → annotate → checkpoints |
| [`4_GSEA.R`](4_GSEA.R) | fgsea on MSigDB Hallmark + the ALS WikiPathways set |
| [`5_Graphs.R`](5_Graphs.R) | PCA, volcanoes, heatmaps, per-gene bar graphs, Venn diagrams |
| [`5b_pathway_QC.R`](5b_pathway_QC.R) | Redundancy check on the Figure 7 pathway panels |

Import is separated from QC so that re-running QC while tuning thresholds does not mint a
new counts checkpoint each time. `1_df_count.R` writes `annotated_counts` and
`vendor_genes`; scripts 2 and 3 both load them and neither writes them back.
`load_checkpoint()` takes the highest version, so the version it reports should match the
"Counts checkpoint" row of `results/qc_report.md`.

Significance cutoffs live only in `0_config.R` (`P_CUT`, `LFC_CUT`); do not redeclare them
per script. QC cutoffs live only in `QC_THRESHOLDS` in `2_QC.R`.

Column contract from script 3 onward is limma's: `logFC`, `AveExpr`, `t`, `P.Value`,
`adj.P.Val`, `B`, plus `ensembl_id` / `entrez_id` / `SYMBOL`. Expression values are voom
log2-CPM (`v$E`), not DESeq2 normalised counts.

## When the data arrives

1. Drop the Plasmidsaurus matrix at `data/raw/WNV_ALS_R01_2026-expression-matrix.tsv`.
2. Run [`1_df_count.R`](1_df_count.R) once.
3. Run [`2_QC.R`](2_QC.R) and read `results/qc_report.md`. Retune `QC_THRESHOLDS` and re-run
   as often as needed — it writes no checkpoint.
4. Fill in `drop_samples <- character(0)` in [`3_DE_limma.R`](3_DE_limma.R) with any samples
   to exclude, then run it.
5. Run scripts 4 and 5.

`data/r_objects/` must hold no checkpoints before the first real run — `load_checkpoint()`
takes the highest version number, so a stale object would be picked up silently.
