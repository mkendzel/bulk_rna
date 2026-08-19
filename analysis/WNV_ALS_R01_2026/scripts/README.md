# WNV_ALS_R01_2026

Human brain organoid bulk RNA-seq: ALS (C9orf72) vs healthy control lines, challenged with
IFN-beta, IFN-gamma, or West Nile virus.

Run everything with the **repo root** (`bulk_rna/`) as the working directory.

## Design

40 samples in two independent experiments.

| Line | Sheet name | Experiment | Mock | IFNb | IFNg | WNV | n |
|---|---|---|---|---|---|---|---|
| `C9` | C9_1 ALS spinal organoids | spinal | 3 | 3 | 3 | 3 | 12 |
| `Ctrl3` | Ctrl_3 healthy control | spinal | 2 | 3 | 3 | 3 | 11 |
| `Ctrl2` | Ctrl_2 healthy control | spinal | 3 | **1** | **1** | 3 | 8 |
| `Cort` | Cortical organoids, uninfected | cort | 3 | 3 | 3 | - | 9 |

The cortical organoids are a separate experiment. They get their own `DGEList`, normalisation
factors and dispersion estimates. Nothing is pooled across the two and spinal and cort logFC /
NES values are not on a shared scale.

Two structural limitations to consider before reading results.

- **Ctrl2 IFNb and IFNg are n = 1.** The cell-means fit still estimates those group means
  (residual variance is pooled across all 12 groups and moderated by `eBayes`), but the six
  contrasts touching them are unreplicated. They carry `min_n = 1` in the contrast registry. 
  Triangles are used in the graphs to showcase this n of 1.

- **The genotype comparison is n = 1 donor.** Replicates are organoid batches from one ALS iPSC
  line against two control lines, so "C9 vs control" has a nonfixable or estimatable donor effect.


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

The vendor matrix does not carry the coded names. Its columns are `FCMSVB_<idx>_cpm` /
`FCMSVB_<idx>_count`, where `<idx>` is the trailing global index of the plasmid code —
`FCMSVB_10` is `C9-TX1_10` is `C9_WNV_1`. `vendor_key()` in [`0_config.R`](0_config.R) reads the
vendor's own idx → code map out of the per-sample filenames inside `data/raw/FCMSVB_results.zip`
(`FCMSVB_mapping-stats/FCMSVB_<idx>_<plasmid>.tsv`), and script 1 asserts it matches `sample_map`
before matching columns in three tiers: vendor key, exact plasmid code, trailing index. Unmapped
columns fail `stopifnot` rather than being dropped.

## Raw files

Only two files in `data/raw/` are used:

| File | Used for |
|---|---|
| `FCMSVB-expression-matrix.tsv` | counts and CPM, 40 samples |
| `FCMSVB_results.zip` | sample-name map, alignment read stats, MultiQC / PCA / correlation reports |

The results zip is read in place with `unzip(list = TRUE)` and `unz()`; nothing is extracted into
the repo. All of `data/raw/` is gitignored.

## Contrasts

25 contrasts, built in `contrast_registry` ([`0_config.R`](0_config.R)). Scripts 4–7 join against
that table on `name` to recover grouping metadata and plot labels.

| type | n | pattern |
|---|---|---|
| `within_line` | 11 | `<line>_<stim> - <line>_Mock`, all four lines |
| `genotype` | 8 | `C9_<stim> - Ctrl3_<stim>` and `C9_<stim> - Ctrl2_<stim>` |
| `interaction` | 6 | `(C9_<stim> - C9_Mock) - (<ref>_<stim> - <ref>_Mock)` |

The interaction family is why this pipeline uses limma-voom rather than DESeq2: `makeContrasts` on
a `~ 0 + condition` cell-means design expresses a difference of differences directly, whereas
`lfcShrink(type = "apeglm")` accepts only a single coefficient. Those six contrasts are also the
weakest test by construction — they sum variance across four group means, giving a median standard
error of ~0.44 against ~0.28 for a within-line contrast.

Not built: a cortical-vs-spinal mock comparison, which would confound tissue with batch.

## Scripts

| | |
|---|---|
| [`0_config.R`](0_config.R) | Paths, cutoffs, `sample_map`, `contrast_registry`, palettes. Sourced by all others. |
| [`1_df_count.R`](1_df_count.R) | Import → rename → checkpoint. Run once per delivery. |
| [`2_QC.R`](2_QC.R) | 14 QC checks → `results/qc_report.md` + `figures/qc/`. Writes no checkpoint. |
| [`3_DE_limma.R`](3_DE_limma.R) | `drop_samples` → split by experiment → voom/lmFit/eBayes → annotate → checkpoints |
| [`4_GSEA.R`](4_GSEA.R) | fgsea on MSigDB Hallmark + the ALS WikiPathways set → checkpoints + TSVs. Draws nothing. |
| [`5_Graphs.R`](5_Graphs.R) | Overview, DE and expression figures: PCA, sample correlation, heatmaps, volcanoes, MA, p-value histograms, DE counts, interaction scatters, Venn/UpSet, gene panels |
| [`6_GSEA_figures.R`](6_GSEA_figures.R) | GSEA figures: NES heatmap, bubble, per-contrast lollipops, running-ES curves, leading-edge heatmaps, custom gene set NES |
| [`7_pathway_QC.R`](7_pathway_QC.R) | Leading-edge redundancy and keyword QC on the significant Hallmark calls, per contrast |
| [`8_composition.R`](8_composition.R) | Marker-module scores for cell composition / maturity → `figures/composition/` + `results/composition_module_stats.tsv` |
| [`9_composition_adjusted.R`](9_composition_adjusted.R) | Whether the IFN result is separable from composition: per-sample ISG regression, control-only calibration, limma with/without the astrocyte covariate |
| [`make_contrast_table.R`](make_contrast_table.R) | Standalone: renders `contrast_registry` to `notes/WNV_ALS_R01_2026_contrasts.pdf` |
| [`make_model_schematic.R`](make_model_schematic.R) | Standalone: documentation figure of the limma-voom model, sized from the checkpoints |

Order of operations:

1. Run [`1_df_count.R`](1_df_count.R) once per delivery. `data/r_objects/` must hold no checkpoints
   before the first run — `load_checkpoint()` takes the highest version number, so a stale object
   would be picked up silently.
2. Run [`2_QC.R`](2_QC.R) and read `results/qc_report.md`. Retune `QC_THRESHOLDS` and re-run as
   often as needed; it writes no checkpoint.
3. Set `drop_samples` in [`3_DE_limma.R`](3_DE_limma.R) to any samples to exclude, then run it.
4. Run scripts 4–9.

Import is separate from QC so that re-running QC while tuning thresholds does not mint a new counts
checkpoint. `1_df_count.R` writes `annotated_counts` and `vendor_genes`; scripts 2 and 3 both load
them and neither writes them back. The version `load_checkpoint()` reports should match the "Counts
checkpoint" row of `results/qc_report.md`.

Figures belong in `5_Graphs.R` and `6_GSEA_figures.R` (QC figures in `2_QC.R` are the exception).
Scripts 3 and 4 compute and checkpoint only, so a plot can be retuned without re-fitting or
re-permuting.

**`fgsea()` is never called outside `4_GSEA.R`.** It sets no seed and its p-values are
permutation-based, so a re-run would not reproduce `results/gsea_hallmark_all.tsv` and the padj
drawn on a figure would contradict the table. Script 6 rebuilds the rank vectors with
`rank_stats()` — asserted at run time to match `run_gsea()`'s recipe for all 25 contrasts — and uses
only `fgsea::plotEnrichment()`, which is deterministic.

## QC

[`2_QC.R`](2_QC.R) runs 14 per-sample checks and writes `results/qc_report.md` plus ten plots in
`figures/qc/`. It filters nothing and writes no checkpoint; `drop_samples` lives in
[`3_DE_limma.R`](3_DE_limma.R).

`qc_status_grid.png` gives every sample against every check at a glance;
`qc_metrics_by_sample.png` is the same data as bars with the cutoffs drawn in. The RNAseqQC figures
beside them describe the run as a whole rather than naming samples.

Checks come in three groups: vendor alignment stats from `FCMSVB_results.zip`, count matrix
statistics (library size, complexity, gene detection, biotype composition), and vst-derived
agreement measures (within-condition correlation, PCA centroid distance, replicate deviation). Plot
categories follow the
[RNAseqQC vignette](https://cran.r-project.org/web/packages/RNAseqQC/vignettes/introduction.html);
`plot_chromosome()` is unused because the vendor matrix carries no gene coordinates.

**Every cutoff lives in the `QC_THRESHOLDS` tribble at the top of `2_QC.R`**, with non-threshold
knobs in `QC_PARAMS` beside it. All values are absolute; nothing else in the script hardcodes a
limit. The vst-derived cutoffs (`within_condition_cor`, `pca_centroid_dist`, `median_abs_M`) are
calibrated against this dataset and will need retuning on a different run. `mapping_rate`,
`dedup_rate`, `protein_coding_pct` and `rRNA_pct` are sentinels set far from anything these samples
produce.

As of the 2026-08 run: 506 pass, 30 warn, 20 fail; seven samples fail at least one check and
`drop_samples` is empty.

## Figures

All figures go through `save_fig()` in [`R/save_fig.R`](../../../R/save_fig.R): PNG at 300 dpi by
default (`FIG_FORMATS`), `formats = FIG_FORMATS_PUB` for the pdf+tiff manuscript panels. Never call
`ggsave()` or `pheatmap(filename = )` directly — `save_fig()` is what sets `bg = "white"` and
handles pheatmap gtables and UpSetR base-device plots.

`figures/` is gitignored. Roughly 170 files, laid out family first, experiment second:

| Directory | Contents |
|---|---|
| `qc/` | `2_QC.R` checks + voom mean-variance from script 3 |
| `overview/` | `pca_<e>`, `pca_pc34_<e>`, `pca_scree_<e>`, `sample_cor_<e>`, `heatmap_top50_variable_<e>`, `heatmap_log2FC_<e>`, `de_counts_<e>` |
| `de/volcano/<e>/` | `volcano_<contrast>` — one per contrast; `de/volcano/volcano_grid_<e>` for the comparative read |
| `de/ma/<e>/` | `ma_<contrast>` — one per contrast; `de/ma/ma_grid_<e>` |
| `de/diagnostics/` | `pvalue_hist_<e>` — faceted, annotated with the enrichment ratio and π₀ |
| `de/overlap/` | `venn_<stim>`, `upset_genotype_<dir>`, `upset_within_line_<dir>` |
| `interaction/` | `interaction_<contrast>` + `interaction_overview` |
| `panels/<e>/` | `panel_<setname>` — faceted multi-gene expression, one file per panel |
| `gsea/hallmark/` | `nes_heatmap_<e>`, `bubble_<e>` |
| `gsea/hallmark/lollipop/<e>/` | `top_pathways_<contrast>` — one per contrast |
| `gsea/hallmark/curves/<e>/` | `curves_<contrast>` — top 5 running-ES curves per contrast, one patchwork |
| `gsea/hallmark/leading_edge/<e>/` | `le_<PATHWAY>` — fixed 4-pathway list, genes × samples |
| `gsea/pathway_qc/<e>/` | `le_jaccard_<contrast>` — leading-edge overlap from script 7 |
| `gsea/als/` | `als_NES_<e>` |
| `composition/` | `composition_modules_<e>`, `composition_markers_<e>`, `composition_adjusted_isg_<e>` |

Every per-experiment figure is built once per entry in `EXPERIMENTS` and saved with the experiment
in the filename.

Every figure is written even when there is nothing to draw — `empty_plot()` in
[`R/plot_guards.R`](../../../R/plot_guards.R) is the placeholder, and nothing in scripts 5–7 uses
`next` or `return(NULL)` in place of a figure. Likewise all 25 contrasts appear on every
per-experiment axis: discrete scales are built from `contrast_labels(e)` with `drop = FALSE`, and
the GSEA bubble draws non-significant cells as grey crosses rather than filtering rows away.
Without those, ggplot's dropping of unused discrete levels and `facet_grid(scales = "free_x")`'s
dropping of empty panels would make a contrast with nothing significant silently disappear — which
the six interaction contrasts would do, since none has a gene past `adj.P.Val < P_CUT`.

## Conventions

- Significance cutoffs live only in `0_config.R` (`P_CUT`, `LFC_CUT`); QC cutoffs live only in
  `QC_THRESHOLDS` in `2_QC.R`. Do not redeclare either per script.
- Column contract from script 3 onward is limma's: `logFC`, `AveExpr`, `t`, `P.Value`,
  `adj.P.Val`, `B`, plus `ensembl_id` / `entrez_id` / `SYMBOL`.
- Expression values are voom log2-CPM (`v$E`), not DESeq2 normalised counts.
