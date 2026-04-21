# Repository for r_objects generated during data processing and analysis of Experiment LSS7T8 using the **DESeq2 Pipeline**

## `plasmid_counts(all)`

**File:** `data/LSS7T8/r_objects/plasmid_counts(all).rds`
### Description
Raw Plasmidsaurus count matrix, 18 samples across three treatment groups (control, CL, DOPC) and two tissues (liver, spleen). Columns renamed from instrument-assigned
IDs (`LSS7T8_1` through `LSS7T8_18`) to `{treatment}_{tissue}_{n}`.
### Sample Layout

| Original ID   | Renamed Label      | Treatment | Tissue | Sample # |
|---------------|--------------------|-----------|--------|----------|
| LSS7T8_1      | control_liver_1    | Control   | Liver  | 1        |
| LSS7T8_2      | control_liver_2    | Control   | Liver  | 2        |
| LSS7T8_3      | control_liver_3    | Control   | Liver  | 3        |
| LSS7T8_4      | control_spleen_1   | Control   | Spleen | 1        |
| LSS7T8_5      | control_spleen_2   | Control   | Spleen | 2        |
| LSS7T8_6      | control_spleen_3   | Control   | Spleen | 3        |
| LSS7T8_7      | cl_liver_4         | CL        | Liver  | 4        |
| LSS7T8_8      | cl_liver_5         | CL        | Liver  | 5        |
| LSS7T8_9      | cl_liver_6         | CL        | Liver  | 6        |
| LSS7T8_10     | cl_spleen_4        | CL        | Spleen | 4        |
| LSS7T8_11     | cl_spleen_5        | CL        | Spleen | 5        |
| LSS7T8_12     | cl_spleen_6        | CL        | Spleen | 6        |
| LSS7T8_13     | dopc_liver_7       | DOPC      | Liver  | 7        |
| LSS7T8_14     | dopc_liver_8       | DOPC      | Liver  | 8        |
| LSS7T8_15     | dopc_liver_9       | DOPC      | Liver  | 9        |
| LSS7T8_16     | dopc_spleen_7      | DOPC      | Spleen | 7        |
| LSS7T8_17     | dopc_spleen_8      | DOPC      | Spleen | 8        |
| LSS7T8_18     | dopc_spleen_9      | DOPC      | Spleen | 9        |
### Processing Notes
- Renaming performed defensively via `intersect()` against existing column names;
  missing columns produce no error and are flagged by `setdiff()` post-hoc.
- Duplicate column name check confirmed `FALSE` after renaming.
- Sample numbering is global (1--9) rather than per-group (1--3), reflecting
  instrument run order.
### QC
Per-sample QC computed across total reads, detected genes (count > 0), and each
sample's percentage of total library depth. No samples were removed. Two samples
(`cl_spleen_6`, `dopc_spleen_9`) were flagged as low-count candidates in a prior
run but retained here; exclusion lines are preserved in the script as
commented-out code.
### Contents (as saved)
Full renamed count matrix, all 18 samples, no rows or columns removed.


---
---
## `sample_info` -- Sample Metadata Table

**File:** `data/LSS7T8/r_objects/sample_info.rds`

### Description
Data frame of sample-level metadata parsed from column names of `counts` via regex. One row per sample, three columns: `sample`, `treatment`, `tissue`.

### Processing Notes
- `treatment` and `tissue` extracted by regex stripping; no lookup table used.
- Both columns converted to factors; `treatment` releveled with `"control"` as reference.
- Row names set to sample names for direct alignment with count matrix columns.
---
## `voom_liver` -- voom-Transformed Counts, Liver

**File:** `data/LSS7T8/r_objects/voom_liver.rds`

### Description
EList object produced by `limma::voom()` on the filtered, TMM-normalized liver DGEList. Contains observation-level precision weights used by `lmFit`.

### Processing Notes
- DGEList built from `counts_liver` with counts rounded to integer.
- Low-expression genes removed via `filterByExpr()` (default thresholds, grouped by `treatment`); library sizes recalculated after filtering (`keep.lib.sizes = FALSE`).
- TMM normalization applied via `calcNormFactors()` before voom transformation.
- Design matrix: `~ treatment`, reference level `"control"`.

---
## `res_voom_liver` -- voom/limma Results, Liver

**File:** `data/LSS7T8/r_objects/res_voom_liver.rds`

### Description
Named list of `topTable` data frames, one per pairwise treatment comparison for liver, from a voom/limma pipeline.

### Comparisons
- `treatmentcl` -- CL vs. control, from `colnames(design)`
- `treatmentdopc` -- DOPC vs. control, from `colnames(design)`
- `treatment_dopc_vs_cl` -- added via `makeContrasts(treatmentdopc - treatmentcl)`

### Processing Notes
- All comparisons extracted with `number = Inf, sort.by = "none"` to return the full gene list in a consistent order.
- The DOPC vs. CL contrast requires a second `contrasts.fit`/`eBayes` call (`fit2`) because it cannot be expressed as a single design matrix coefficient.
- Note that coefficient names for the reference-level comparisons follow limma's convention (`treatmentcl`, `treatmentdopc`)


---
## `res_voom_liver_annotated` -- voom/limma Results with Gene Annotation, Liver

**File:** `data/LSS7T8/r_objects/res_voom_liver_annotated.rds`

### Description
Annotated version of `res_voom_liver`. Each comparison is a tibble with limma results joined to gene symbols and Entrez IDs from `mouse_gene_map`.

### Processing Notes
- ENSEMBL version suffixes stripped via `sub("\\..*$", "", ...)` before joining.
- Many-to-one mappings resolved by `distinct(ENSEMBL, .keep_all = TRUE)` on `mouse_gene_map` prior to join.
- Genes with no annotation match retain `NA` for `SYMBOL` and `ENTREZID`.