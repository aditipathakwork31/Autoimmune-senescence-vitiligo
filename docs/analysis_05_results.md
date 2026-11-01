# Analysis 05 results (pre-specified in analysis_plan_05.md)

Run date: 2026-11-01. No rule, marker list or threshold was changed after seeing results.
The live GEO lookup failed, so the sample titles embedded in the script were used.

## Data
- GSE298871: 17 whole-skin 10x samples (6 healthy, 6 lesional, 5 non-lesional); 24,746 cells, 22,462 after filtering; 9.2% unassigned.
- Cells per class: keratinocyte 6,962; endothelial 5,279; fibroblast 3,519; mural 2,534; T cell 1,295; myeloid 479; melanocyte 330; mast 1.
- Descriptive check (not pre-specified): melanocytes were nearly absent from lesional samples (3 cells in 6 samples), present in non-lesional (243) and healthy (84).

## Part 1: whole-sample replication
- S0 (SenMayo): slope 0.58 (95% CI 0.12 to 1.05), p = 0.018: coupled.
- S1 (arrest panel): slope 0.50 (-0.01 to 1.01), p = 0.055: inconclusive (not clearly coupled, not clearly weaker).
- Outcome: P1-inconclusive. Paired subset (5 patients): intervals too wide to interpret.

## Part 2: cell types
- Not testable (fewer than 8 eligible samples): T cell (6), myeloid (2). Fibroblast paired comparison not testable (3 pairs).
- Keratinocyte (16 samples): S0 slope 0.66 (0.18 to 1.15), p = 0.012, BH 0.047: coupled. S1 0.35 (-0.23 to 0.94), p = 0.21: inconclusive.
- Fibroblast (14 samples): S0 0.53 (-0.13 to 1.18), p = 0.105: inconclusive. S1 0.14 (-0.54 to 0.82), p = 0.65: not coupled.
- Flags: F1 not testable; F2 yes; F3 no. Interpretation (fixed in the plan): F2 yes, the coupling is not confined to immune cells.

## Caveats
- The immune versus structural comparison could not be made because T cells and myeloid cells were too few.
- F2 rests on one cell type at BH-adjusted p = 0.047; the S1 keratinocyte result is inconclusive by a narrow margin.
- SenMayo is largely inflammatory secreted factors: coupling in keratinocytes does not show senescence.
- The analysis 04 conclusion about the arrest panel is weakened: here its slope is close to SenMayo's.
- Samples from the same patient are treated as independent in the unpaired models; unit size (50 to about 1,000 cells) was not examined.
- GSE203262 (Part 3) was not analysed: its lesional/non-lesional barcode mapping is undocumented.

Full tables: results/05_singlecell_vitiligo/
