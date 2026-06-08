# Repository for r_objects generated during data processing and analysis

## Overview
There are two experiments whose processed data is stored in this repository: R2SDHF and C9orf72
- R2SDHF - Shilu Malakar's bulk RNA-seq data (sequenced and aligned via Plasmidsaurus)
- C9orf72 - Jimena Andersen's bulk RNA-seq data (unknown sequence or alignment company).

----
## **Experiment: R2SDHF**
### `R2SDHF.rds`
**Source counts:** `counts_R2SDHF`
**Metadata:** `meta_R2SDHF` (grouping variable: `condition`)

**Pipeline:** edgeR/limma-voom

- Built `DGEList`, filtered low-expression genes via `filterByExpr` (group = `condition`), recalculated library sizes.
- TMM normalization (`calcNormFactors`).
- Design: `~ condition` from `meta_R2SDHF`.
- `voom` → `lmFit` → `eBayes`.

**Gene mapping:** Converted Ensembl IDs to HGNC symbols via `org.Hs.eg.db`. Dropped unmapped (`NA`) and duplicate symbols.

**Saved object:** `v_R2SDHF_sym` — the voom-normalized log-CPM expression matrix (`v_R2SDHF$E`) subset to uniquely mapped genes, with gene symbols as rownames. This is the expression matrix only, not the full `EList`.

**Path:** `data/Anderson-Suthar_collab/Als_org_Jimena/R2SDHF.rds`

### R2SDHF objects (revised pipeline)

**Source counts:** `counts_R2SDHF`
**Metadata:** `meta_R2SDHF` (grouping variable: `Group`)

**Pipeline:** edgeR/limma-voom

- Design: `~ 0 + Group` (no intercept, group-means parameterization). Column names set to factor levels of `Group`.
- Built `DGEList`, filtered low-expression genes via `filterByExpr` (using design matrix), recalculated library sizes.
- TMM normalization (`calcNormFactors`).
- `voom` → `lmFit`.

**Contrasts:**
- `TX12_vs_Mock` = TX_12 − Mock_0
- `TX24_vs_Mock` = TX_24 − Mock_0
- `TX48_vs_Mock` = TX_48 − Mock_0
- `TX48_vs_TX12` = TX_48 − TX_12

Applied via `contrasts.fit` → `eBayes`.

**Note:** Gene IDs remain as-is (Ensembl); no symbol mapping in this pipeline. The earlier `R2SDHF.rds` (expression matrix with HGNC symbols) was from a prior version using `~ condition` with an intercept design.

#### Saved objects

| File | Object | Description |
|------|--------|-------------|
| `Als_org_Jimena/fit2(R2SDHF).rds` | `fit2` | eBayes fit with all four contrasts |
| `Als_org_Jimena/voom(R2SDHF).rds` | `v` | Full voom `EList` (Ensembl IDs) |
| `r_objects/res_TX12(R2SDHF).rds` | `res_TX12` | `topTable` — TX_12 vs Mock_0 |
| `r_objects/res_TX24(R2SDHF).rds` | `res_TX24` | `topTable` — TX_24 vs Mock_0 |
| `r_objects/res_TX48(R2SDHF).rds` | `res_TX48` | `topTable` — TX_48 vs Mock_0 |
| `r_objects/res_TX48_vs_TX12(R2SDHF).rds` | `res_TX48_vs_TX12` | `topTable` — TX_48 vs TX_12 |

### R2SDHF contrasts (fgsea)

**Gene sets:** MSigDB Hallmark collection (v2025.1, human, symbols) — `h.all.v2025.1.Hs.symbols.gmt`
**Input:** `topTable` results from revised R2SDHF pipeline (gene symbols as rownames)

**Method:** `fgsea` with 10,000 permutations. Ranking statistic: moderated t-statistic (`tt$t`). Duplicate gene names resolved by keeping the entry with highest `AveExpr` before ranking.

**Assumes** gene symbol mapping was applied to `topTable` rownames upstream (not shown in the revised pipeline code — see note on prior entry).

#### Saved objects

| File | Object | Contrast |
|------|--------|----------|
| `r_objects/fgsea_TX12_vs_Mock(R2SDHF).rds` | `fgsea_TX12` | TX_12 vs Mock_0 |
| `r_objects/fgsea_TX24_vs_Mock(R2SDHF).rds` | `fgsea_TX24` | TX_24 vs Mock_0 |
| `r_objects/fgsea_TX48_vs_Mock(R2SDHF).rds` | `fgsea_TX48` | TX_48 vs Mock_0 |
| `r_objects/fgsea_TX48_vs_TX12(R2SDHF).rds` | `fgsea_TX48_vs_TX12` | TX_48 vs TX_12 |



----
## **Experiment: C9ors72 current pipeline (hSpS, Ctrx, and batch)**

### Counts & metadata setup (hSpS, Ctrx, and batch)

**Source files:**
- Expression: `Als_org_Jimena/102025_lines1-6,c9_vs_ctrl,d30,50,75,120_allfeature_counts.xlsx`
- Metadata: `Als_org_Jimena/BulkSeqTargetList_tosend.xlsx`

**Expression matrix cleanup:** Set `Geneid` as rownames. Stripped leading `./` and trailing `_S*` from column names to get clean sample IDs.

**Metadata subset:** Filtered to hSpS organoids on Ctrx coating, ALS vs control lines only. Retained columns: `Admera_Health_ID`, `Patient`, `AgeOrg`, `Line`, `Replicate`, `Sex`, `Batch` (renamed from `d0`).

**Factor setup:**
- `Line`: control (ref) vs ALS
- `AgeOrg`, `Patient`, `Sex`, `Batch`: all factored

**Derived variable — `Stage`:** Binned `AgeOrg` into developmental stages:
- d30: ages 30, 31, 34, 35
- d50: ages 50, 52
- d75: ages 75, 76, 77
- d120: age 120
- Ordered factor: d30 → d50 → d75 → d120

**Counts matrix:** Subset of `expr` columns matching `meta_filtered$Admera_Health_ID`.

**QC:** Per-sample total reads, detected genes, percent of total library. Per-sample gene-level summary stats. Not saved — exploratory only.

#### Saved objects

| File | Object | Description |
|------|--------|-------------|
| `Als_org_Jimena/counts(hSpS_Ctrx_batch).rds` | `counts` | Raw count matrix, hSpS/Ctrx samples only |
| `Als_org_Jimena/meta_filtered(hSpS_Ctrx_batch).rds` | `meta_filtered` | Metadata with factors and `Stage` variable |

### DE analysis: ALS vs Control by Stage (hSpS, Ctrx, and batch)

**Input:** `counts` and `meta_filtered` from hSpS Ctrx batch setup.

**Pipeline:** edgeR/limma-voom with repeated measures

- Filtered low-expression genes via `filterByExpr` (group = `Line`), recalculated library sizes.
- TMM normalization.
- Design: `~ Sex + Batch + Stage * Line` → adjusts for sex and batch; tests ALS vs control effect at each stage via interaction terms. Reference levels: control, d30.
- `duplicateCorrelation` run twice (recommended iteration) with blocking on `Patient` to account for repeated measures across stages within the same patient.
- `voom` → `lmFit` (both with patient blocking and consensus correlation) → `eBayes`.

**Contrasts:** Extract the ALS vs control effect at each stage from the interaction model:
- `d30_ALS` = `LineALS` (main effect at reference stage)
- `d50_ALS` = `LineALS` + `Staged50.LineALS`
- `d75_ALS` = `LineALS` + `Staged75.LineALS`
- `d120_ALS` = `LineALS` + `Staged120.LineALS`

#### Saved objects

| File | Object | Description |
|------|--------|-------------|
| `Als_org_Jimena/fit2(hSpS_Ctrx_batch).rds` | `fit2` | eBayes fit with all four stage-specific ALS contrasts |
| `Als_org_Jimena/voom(hSpS_Ctrx_batch).rds` | `v` | Full voom `EList` (with patient correlation) |
| `r_objects/res_d30_ALS_vs_Ctrl(Batch).rds` | `res_d30` | `topTable` — ALS vs Ctrl at d30 |
| `r_objects/res_d50_ALS_vs_Ctrl(Batch).rds` | `res_d50` | `topTable` — ALS vs Ctrl at d50 |
| `r_objects/res_d75_ALS_vs_Ctrl(Batch).rds` | `res_d75` | `topTable` — ALS vs Ctrl at d75 |
| `r_objects/res_d120_ALS_vs_Ctrl(Batch).rds` | `res_d120` | `topTable` — ALS vs Ctrl at d120 |

### fgsea results: ALS vs Control by Stage (hSpS, Ctrx, and batch)

**Gene sets:** MSigDB Hallmark collection (v2025.1, human, symbols) — `h.all.v2025.1.Hs.symbols.gmt`
**Input:** `topTable` results from hSpS Ctrx batch DE analysis (stage-specific ALS vs Ctrl contrasts)

**Method:** Same `run_gsea` helper as R2SDHF analysis. `fgsea` with 10,000 permutations. Ranking statistic: moderated t-statistic. Duplicates resolved by keeping highest `AveExpr`.

#### Saved objects

| File | Object | Contrast |
|------|--------|----------|
| `r_objects/fgsea_d30_ALS_vs_Ctrl(Batch).rds` | `gsea_d30` | ALS vs Ctrl at d30 |
| `r_objects/fgsea_d50_ALS_vs_Ctrl(Batch).rds` | `gsea_d50` | ALS vs Ctrl at d50 |
| `r_objects/fgsea_d75_ALS_vs_Ctrl(Batch).rds` | `gsea_d75` | ALS vs Ctrl at d75 |
| `r_objects/fgsea_d120_ALS_vs_Ctrl(Batch).rds` | `gsea_d120` | ALS vs Ctrl at d120 |

### clusterProfiler GSEA results: ALS vs Control by Stage (hSpS Ctrx batch)

**Gene sets:** MSigDB Hallmark collection (v2025.1, human, symbols) — same GMT as fgsea runs.
**Input:** Same `topTable` results from hSpS Ctrx batch DE analysis.

**Method:** `clusterProfiler::GSEA`. GMT converted to `TERM2GENE` (t2g) data frame format. Ranking statistic: moderated t-statistic, deduplicated (first occurrence kept after sort). `pvalueCutoff = 1` (retain all pathways), `eps = 0` (exact p-values).

**Why both fgsea and clusterProfiler:** fgsea objects are plain data frames; clusterProfiler returns an S4 `gseaResult` object compatible with `enrichplot` visualization functions (`dotplot`, `ridgeplot`, `gseaplot2`, etc.).

#### Saved objects

| File                                            | Object    | Contrast            |
| ----------------------------------------------- | --------- | ------------------- |
| `r_objects/cp_gsea_d30_ALS_vs_Ctrl(Batch).rds`  | `cp_d30`  | ALS vs Ctrl at d30  |
| `r_objects/cp_gsea_d50_ALS_vs_Ctrl(Batch).rds`  | `cp_d50`  | ALS vs Ctrl at d50  |
| `r_objects/cp_gsea_d75_ALS_vs_Ctrl(Batch).rds`  | `cp_d75`  | ALS vs Ctrl at d75  |
| `r_objects/cp_gsea_d120_ALS_vs_Ctrl(Batch).rds` | `cp_d120` | ALS vs Ctrl at d120 |


### DE & GSEA pipeline (OUTDATED) (no batch correction)

**Superseded by:** hSpS Ctrx batch pipeline. This version lacks `Batch` in the design matrix.

**Source data:**
- Counts: `Als_org_Jimena/counts(hSpS_Ctrx).rds`
- Metadata: `Als_org_Jimena/meta_filtered(hSpS_Ctrx).rds`

**Pipeline:** edgeR/limma-voom with repeated measures

- Filtered low-expression genes via `filterByExpr`, recalculated library sizes. TMM normalization.
- Design: `~ Sex + Stage * Line` — adjusts for sex; tests ALS vs control at each stage via interaction. Missing `Batch` covariate (corrected in later pipeline).
- `duplicateCorrelation` run twice with blocking on `Patient`.
- `voom` → `lmFit` (both with patient blocking and consensus correlation) → `eBayes`.

**Contrasts:** ALS vs Control at each stage:
- d30: `LineALS`
- d50: `LineALS + Staged50.LineALS`
- d75: `LineALS + Staged75.LineALS`
- d120: `LineALS + Staged120.LineALS`

#### DE results

`topTable` data frames (one row per gene): `logFC`, `AveExpr`, `t`, `P.Value`, `adj.P.Val`, `B`. Positive logFC = higher in ALS.

| File | Contrast |
|------|----------|
| `r_objects/res_d30_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d30 |
| `r_objects/res_d50_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d50 |
| `r_objects/res_d75_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d75 |
| `r_objects/res_d120_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d120 |

#### fgsea results

MSigDB Hallmark (v2025.1, human, symbols). Ranked by moderated t-statistic, 10,000 permutations. Duplicates resolved by highest `AveExpr`.

| File | Contrast |
|------|----------|
| `r_objects/fgsea_d30_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d30 |
| `r_objects/fgsea_d50_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d50 |
| `r_objects/fgsea_d75_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d75 |
| `r_objects/fgsea_d120_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d120 |

#### clusterProfiler GSEA results

Same gene sets and ranking as fgsea. `gseaResult` S4 objects for use with `enrichplot` (`dotplot`, `ridgeplot`, `gseaplot2`, `emapplot`, `cnetplot`, etc.). `pvalueCutoff = 1`, `eps = 0`.

| File | Contrast |
|------|----------|
| `r_objects/cp_gsea_d30_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d30 |
| `r_objects/cp_gsea_d50_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d50 |
| `r_objects/cp_gsea_d75_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d75 |
| `r_objects/cp_gsea_d120_ALS_vs_Ctrl.rds` | ALS vs Ctrl at d120 |
