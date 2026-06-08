# Spinal Cord Organoid Datasets

Datasets were generated from spinal cord organoids.

## Timepoints

- **D30:** Postmitotic motor neurons, some other young neurons and progenitors
- **D50:** Postmitotic neurons, some progenitors
- **D75:** Postmitotic neurons, starting generation of astrocytes
- **D120:** Postmitotic neurons, astrocytes, oligodendrocyte lineage cells

## hiPSC Lines

Lines are derived from:

- Healthy controls
- C9orf72 expansion repeat patients
- SOD1 A5V mutation patients

---

## Bulk Dataset

DE lists show C9orf72 vs healthy control:

- Negative fold change = down in C9orf72
- Positive fold change = up in C9orf72
- Comparison of 6 healthy controls and 6 C9orf72 patients (at least 2 batches per line)
- Separate files for D30, D50, and D75

---

## Single Cell Datasets

### Neurons

- **Timepoint:** D50
- **Enrichment:** Immunopanned using THY1 and CD24
- **Lines:** 4 healthy control, 4 C9orf72
- Separated into neuronal subtypes (including motor neurons and different interneuron types)

### Astrocytes

- **Timepoint:** D120
- **Enrichment:** Immunopanned using HEPACAM1
- **Lines:** 3 healthy control, 2 C9orf72, 2 SOD1
- **Conditions:**
  - Baseline: no treatment
  - Cytokine treatment: TNFa, IL1a, and C1q for 30 days prior to collection to mimic neuroinflammation
- Not separated into subtypes

### Oligodendrocyte Lineage

- **Timepoint:** D120
- **Enrichment:** Immunopanned using O4
- **Lines:** 3 healthy control, 2 C9orf72, 2 SOD1
- Includes the whole oligodendrocyte lineage (OPCs and mature oligos)
- **Conditions:**
  - Baseline: no treatment
  - Cytokine treatment: TNFa, IL1a, and C1q for 30 days prior to collection to mimic neuroinflammation