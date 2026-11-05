# Analysis 04 results (pre-specified in analysis_plan_04.md)

Run date: 2026-10-28. No rule, gene list or threshold was changed after seeing results.

## Step A: validation in senescent fibroblasts (GSE130100, GSE63577)
- Validated in both datasets: S0 (SenMayo) and S1 (arrest panel: CDKN2A, CDKN2B, GADD45A up; LMNB1 down).
- Dropped: S2 (TANG p53 targets; wrong direction in GSE130100) and S3 (CHICAS RB1 targets;
  significantly lower in senescent cells in GSE63577, so it behaves as a set that falls in senescence).
- CDKN1A was removed from the panel by the plan's rule (it is in the IFN-gamma Hallmark set).

## Step B: coupling to IFN_c (standardised slope)
- S0 (SenMayo): coupled in GSE75819 (0.52), GSE65127 (0.48) and GSE53146 (1.17).
- S1 (panel): not coupled in all three (0.19, 0.16, -1.04; 95% CIs include 0 in every dataset).
- Rule R1 triggered: SenMayo coupled, the validated non-secreted score not coupled.

## Caveats
- Only one non-secreted score (4 genes) was validated; S1 may be driven mainly by LMNB1.
- In GSE75819 and GSE65127 the S1 confidence intervals include SenMayo's slope, so the data do not
  exclude equal coupling.
- SenMayo coupling remains after adjusting for an E2F proliferation score in all three datasets.
- Paired-by-patient SenMayo slopes differ between datasets (GSE65127 0.56, p = 0.0009; GSE75819 0.30, p = 0.30).
- Cohorts are small and bulk tissue; validation was done in fibroblasts only.
- S3's direction was misjudged in the design; it is left dropped, not flipped.

Full tables: results/04_nonsecreted_senescence/

## Update, 2026-11-01 (after analysis 05)

Analysis 05 (docs/analysis_05_results.md) weakens one conclusion drawn above. In the independent
single-cell dataset GSE298871 the arrest panel S1 was inconclusive in the whole-sample test
(slope 0.50, 95% CI -0.01 to 1.01, p = 0.055; SenMayo 0.58, p = 0.018). Rule R1 was triggered here as
pre-specified and that verdict is not changed, but the statement that the coupling reflects inflammatory
secreted content should be read as unresolved: the arrest panel's lack of coupling in analysis 04 was
not replicated as a clear null.
