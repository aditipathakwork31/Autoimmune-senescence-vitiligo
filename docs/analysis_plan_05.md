# Analysis plan 05: which skin cell types carry the interferon-senescence coupling?

Status: FINAL. Written on 29 October 2026 and approved by the author the same day, before any
cell-level data were scored. Once committed, every rule below is fixed. A change after that must be added to the amendment log
at the end, with a date and a reason, and must be made before the affected analysis is run.

## 1. Background and question

Analysis 04 (docs/analysis_plan_04.md, docs/analysis_04_results.md) used bulk skin data. A SenMayo
senescence score (S0) tracked an IFN-γ response score (IFN_c) in all three vitiligo datasets, while a
validated 4-gene arrest panel (S1) did not clearly do so. Bulk tissue cannot show which cells carry
these signals. Question: in single-cell skin data, which cell types show the S0 and S1 coupling to
IFN_c, and does a fourth, independent dataset agree with analysis 04?

## 2. Data

**Primary: GSE298871** (17 whole-skin 10x samples, Cell Ranger 7.2, Miltenyi whole-skin dissociation kit).
- 6 healthy volunteers, 6 lesional (LS), 5 non-lesional (NL) samples.
- Patients with both LS and NL: 486, DER114, DER115, DER119, DER120 (5 pairs). BIO023 has LS only.
- Cells per sample (barcode counts): 465: 542; 486_LS: 1345; 486_NL: 200; 632: 1106; 658: 593;
  BIO023_LS: 1249; CON161: 975; CON171: 1081; CON194: 1690; DER114_LS: 2701; DER114_NL: 2022;
  DER115_LS: 3485; DER115_NL: 2483; DER119_LS: 2699; DER119_NL: 1073; DER120_LS: 505; DER120_NL: 997.
- Raw integer counts. The authors report inflammatory dysregulation in non-lesional skin relative to
  healthy skin, so non-lesional skin is NOT treated as a normal control; healthy samples are a separate group.

**Secondary, conditional: GSE203262** (6 patients, suction-blister cells, Cell Ranger v2.1.0 aggregation).
Each patient archive holds one matrix with barcode suffixes -1 and -2. The GEO text does not state
which suffix is lesional and which is non-lesional. This dataset is used ONLY if the mapping is
documented in the paper's methods or a supplementary table found before any scoring. Labels are never
inferred from expression, interferon signal, or cell counts. If the mapping cannot be documented,
the dataset is excluded and this is reported.

## 3. Step 0: data processing and cell-type assignment (no scores computed in this step)

**Cell filter.** Keep cells with at least 200 detected genes and mitochondrial read fraction
(symbols starting "MT-") below 20%.

**Gene symbols.** Features are matched by symbol; counts of rows sharing a symbol are summed.

**Cell-type assignment by fixed marker lists** (no clustering). For each cell, normalise counts to
10,000 per cell, take log1p, and compute the mean expression of each class's markers. Assign the cell
to the class with the highest mean if that mean is at least 0.5 and exceeds the second-highest class
mean by at least 0.25; otherwise the cell is "unassigned". Markers exclude interferon-response genes
so that assignment does not depend on interferon activity.
- keratinocyte: KRT5, KRT14, KRT15, KRT1, KRT10, KRT17, KRTDAP, SFN
- T cell: CD3D, CD3E, CD2, TRAC, TRBC2, LCK
- myeloid: LYZ, CD68, AIF1, CD14, FCER1G, TYROBP
- fibroblast: COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA
- endothelial: PECAM1, VWF, CLDN5, EGFL7, PLVAP, ACKR1
- melanocyte: PMEL, MLANA, TYR, DCT, MITF, TYRP1
- mast cell: TPSAB1, TPSB2, CPA3, MS4A2
- mural cell: ACTA2, TAGLN, RGS5, MYL9

**Composition check (the only permitted look before scoring).** Report cells per class per sample and
the fraction unassigned. If more than 40% of cells overall are unassigned, the thresholds may be revised
once; the revision is logged in the amendment log before any score is computed.

**Pseudobulk.** For each sample, and for each class, sum counts over its cells. A (sample, class) unit is
"eligible" only if it has at least 50 cells. Eligible units are converted to log2(CPM + 1) after dropping
genes with CPM <= 1 in fewer than 3 units of that class (the same processing as analysis 04).
A "whole-sample" pseudobulk (all filtered cells of a sample) is also built for Part 1.

## 4. Scores (frozen from analysis 04; no change)

Same definitions as scripts/04_nonsecreted_senescence.R:
- IFN_c: Hallmark IFN-γ response minus genes found in SenMayo, Fridman senescence (up) or Reactome SASP.
- S0: SenMayo minus genes in the IFN-γ Hallmark set, scored with GSVA (Gaussian kernel).
- S1: arrest panel (CDKN2A, CDKN2B, GADD45A up; LMNB1 down), mean of z-scored genes, down-gene sign flipped.
  S1 is scored for a class only if at least 3 of its 4 genes have CPM > 1 in that class; otherwise "not scorable".
- E2F: Hallmark E2F targets, GSVA (proliferation control for sensitivity analyses).
- S2 and S3 were dropped in analysis 04 and are not used.
A GSVA score for a class is computed only if at least 8 eligible samples exist for that class and at
least 10 genes of the set are present.

## 5. Part 1: replication in GSE298871 at whole-sample level

Whole-sample pseudobulk, all 17 samples. Models, with IFN_c and each score z-scored over the 17 samples:
- primary: z(score) ~ z(IFN_c) + group (healthy, non-lesional, lesional)
- secondary: the 10 paired LS/NL samples of the 5 paired patients: z(score) ~ z(IFN_c) + status + patient

Classification, as in analysis 04, using the S0 slope in the same model as the reference:
- coupled: slope > 0 and p < 0.05
- not coupled: p >= 0.05 and slope < 0.5 × the S0 slope
- inconclusive: anything else

Outcome (one of three):
- **P1-consistent:** S0 coupled and S1 not coupled.
- **P1-contradicts:** S1 coupled.
- **P1-inconclusive:** anything else.
This outcome is reported alongside, and does not change, the analysis 04 verdict.

## 6. Part 2: cell-type resolution in GSE298871

**Primary cell types:** keratinocyte, T cell, myeloid, fibroblast. All other classes are described
(cell numbers, eligibility) but not tested, with melanocytes noted as expected to be below the minimum.

**Eligibility for a test.**
- coupling test: at least 8 eligible samples for the class;
- paired lesional vs non-lesional test: at least 4 patients with both samples eligible.
A class that does not meet a requirement is reported as "not testable", never tested with less.

**Coupling test, per eligible class:** z(score) ~ z(IFN_c) + group, over the eligible samples, with scores
and IFN_c z-scored within the class. Classification as in Part 1, with the S0 slope in the same class
as the reference.

**Flags (several can apply; each is reported as yes, no or not testable):**
- **F1, coupling in immune cells:** S0 coupled in T cells or myeloid cells.
- **F2, coupling in structural cells:** S0 coupled in keratinocytes or fibroblasts.
- **F3, arrest panel coupled:** S1 coupled in any eligible primary class. If flagged, the proliferation
  sensitivity (E2F added as covariate) is reported next to it.
Interpretation text, fixed in advance:
- F1 yes and F2 no: consistent with the coupling reflecting inflammatory-cell content.
- F2 yes: the coupling is not confined to immune cells.
- F3 yes: notable, and stated with the caveat that arrest genes can follow proliferation or stress changes.
- Everything else: inconclusive.

**Secondary analyses (no verdict).** Per eligible class, lesional vs non-lesional differences in IFN_c, S0 and S1
(paired by patient, at least 4 pairs), and healthy vs non-lesional differences (unpaired, run only with at least 3 eligible samples in each
group, otherwise "not testable"); the E2F-adjusted slopes; per-class eligibility tables.

## 7. Part 3: GSE203262, conditional

Run only if the lesional/non-lesional mapping of the -1 and -2 barcodes is documented before scoring.
Same cell filter, same marker assignment and same eligibility rules (cells are epidermal, so fibroblasts
and many immune classes are expected to be absent). With 6 paired patients, results are descriptive.
Primary classes limited to keratinocyte, T cell and melanocyte. No verdict.

## 8. Multiplicity, power and caveats

- Part 2 reports unadjusted and Benjamini-Hochberg-adjusted p-values across all primary tests;
  classification uses the unadjusted p, as in analysis 04.
- 17 samples (5 pairs) give limited power; "not coupled" is weak evidence of absence, so confidence
  intervals are shown next to every classification.
- Whole-skin dissociation does not sample cell types in proportion to their tissue abundance.
- Sparse single-cell counts give noisy pseudobulk scores for low-abundance genes such as CDKN2A.
- Marker-based assignment is coarse; doublets and ambient RNA are not removed.
- Non-lesional vitiligo skin is inflamed relative to healthy skin, so contrasts with it are conservative.
- Cells from the same patient are not independent; models use sample- or patient-level units only.
- No cell-level p-values are computed.

## 9. What will not be done

- No change to scripts 01-04 or their results.
- No change to the scores, thresholds or rules above after a result is seen.
- No label for any sample or barcode is chosen based on expression values.
- All runs save console output, tables and session information under results/.

## 10. Amendment log

(none yet)

Implementation note, 2026-10-29, before commit: script 05 runs Parts 1 and 2 only. Part 3 (GSE203262) is not
implemented because the barcode-suffix mapping is undocumented.
