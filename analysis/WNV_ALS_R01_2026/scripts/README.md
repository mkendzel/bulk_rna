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

Vendor alignment stats are joined into `qc_summary` as `uniquely_mapped` and drawn as a
third row of `figures/qc/library_size_and_detection.png` with a dashed 5M-read guide.
As of the 2026-08 run, `Cort_Mock_1` (1.6M) and `Ctrl3_IFNb_1` (1.7M) sit far below the
~16M typical, with `Cort_Mock_3` (4.2M) and `Cort_IFNg_2` (7.0M) next. `drop_samples` in
script 1 is empty; set it there after reading the plot.

## Contrasts

25 contrasts, built in `contrast_registry` ([`0_config.R`](0_config.R)). Scripts 2 and 3
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
| [`0_config.R`](0_config.R) | Paths, cutoffs, `sample_map`, `contrast_registry`, palettes. Sourced by 1–3. |
| [`1_QC_DE_limma.R`](1_QC_DE_limma.R) | Import → rename → QC → split by experiment → voom/lmFit/eBayes → annotate → checkpoints |
| [`2_GSEA.R`](2_GSEA.R) | fgsea on MSigDB Hallmark + the ALS WikiPathways set |
| [`3_Graphs.R`](3_Graphs.R) | PCA, volcanoes, heatmaps, per-gene bar graphs, Venn diagrams |

Significance cutoffs live only in `0_config.R` (`P_CUT`, `LFC_CUT`); do not redeclare them
per script.

Column contract from script 1 onward is limma's: `logFC`, `AveExpr`, `t`, `P.Value`,
`adj.P.Val`, `B`, plus `ensembl_id` / `entrez_id` / `SYMBOL`. Expression values are voom
log2-CPM (`v$E`), not DESeq2 normalised counts.

## When the data arrives

1. Drop the Plasmidsaurus matrix at `data/raw/WNV_ALS_R01_2026-expression-matrix.tsv`.
2. Run [`1_QC_DE_limma.R`](1_QC_DE_limma.R) as far as the QC block; read `qc_summary` and
   `figures/qc/library_size_and_detection.png`.
3. Fill in `drop_samples <- character(0)` with any samples to exclude, then run the rest.
4. Run scripts 2 and 3.

`data/r_objects/` must hold no checkpoints before the first real run — `load_checkpoint()`
takes the highest version number, so a stale object would be picked up silently.
