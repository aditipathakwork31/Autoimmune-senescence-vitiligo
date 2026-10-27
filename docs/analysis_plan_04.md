# Analysis plan 04: is the interferon-senescence coupling explained by inflammatory secreted factors?

Status: draft, written on 27 October 2026, before any analysis in this plan has been run.
Everything below is fixed in advance. Any later change must be committed with a date and a reason.

## 1. Question

In scripts 01-03, a SenMayo senescence score tracked an IFN-γ response score in all three vitiligo
datasets, and also in healthy skin. SenMayo contains many inflammatory secreted factors. Does a
senescence score built from cell-cycle-arrest and structural markers, which does not rely on
secreted factors, show the same coupling to interferon?

## 2. Scores (fixed)

All scores are computed per dataset. Genes shared with the IFN-γ Hallmark set (the IFN_c cleaning
rule from scripts 01-03) are removed from every senescence score, and the removed genes are reported.
A gene absent from an array is dropped and reported. A score needs at least 3 present genes
(or the set needs at least 10 present genes for GSVA-scored sets), otherwise it is not computed.

- **S0, SenMayo (reference):** MSigDB SAUL_SEN_MAYO, GSVA, as in scripts 01-03.
- **S1, arrest panel:** up genes CDKN1A, CDKN2A, CDKN2B, GADD45A; down gene LMNB1.
  Score = mean across genes of the within-dataset z-scored expression, with the sign of LMNB1 flipped.
  (TP53 was considered and dropped: its mRNA changes little in senescence.)
- **S2, p53-target senescence:** MSigDB TANG_SENESCENCE_TP53_TARGETS_UP minus TANG_SENESCENCE_TP53_TARGETS_DN,
  each scored with GSVA.
- **S3, RB1-pathway senescence:** MSigDB CHICAS_RB1_TARGETS_SENESCENT, GSVA.
- **Not used:** REACTOME_CELLULAR_SENESCENCE and other sets that contain the SASP sub-pathway,
  because they reintroduce the secreted factors being excluded.
- **Proliferation control (not a senescence score):** MSigDB HALLMARK_E2F_TARGETS, GSVA.

## 3. Step A: validation in senescent fibroblasts (before any vitiligo analysis)

Purpose: show that each score can detect senescence at all. Without this, a null result in
vitiligo cannot be interpreted.

- Datasets: at least two public human-fibroblast datasets (single-channel microarray or RNA-seq)
  with separate proliferating and senescent groups, from different inducers or different studies,
  at least 3 samples per group. Two-colour ratio arrays are excluded. Candidates, to be confirmed
  (sample counts and availability of expression or count tables checked first): GSE130100
  (RAS-induced senescence, BJ-TERT, RNA-seq, 3 control and 3 senescent samples). The second dataset
  is GSE63577 (replicative senescence, RNA-seq, five fibroblast lines); see the amendment log of
  2026-10-28 for how it is used.
- RNA-seq processing: raw counts, genes with low counts filtered, converted to log2 CPM, then
  scored exactly like microarray data (GSVA Gaussian kernel for sets; z-scores for the arrest panel).
- Test per dataset: senescent versus proliferating difference in each score (limma or t-test).
- A score is **validated** if senescent > proliferating (positive direction, p < 0.05) in at least
  2 of the validation datasets.
- **If SenMayo (S0) is not validated, stop.** The validation method is then unreliable and the
  plan is revisited before anything else.
- A score that is not validated is dropped from Step B and reported as dropped.
- Gene-list changes are allowed only here, only for (a) absence from the array, or (b) a gene
  moving in the opposite direction in at least 2 validation datasets. Each change is logged in git.
- **After Step A the gene lists are frozen.** No change after any vitiligo result is seen.

## 4. Step B: coupling to interferon in the vitiligo datasets

Same datasets and models as scripts 01-03:

- GSE75819: score ~ IFN_c + batch + status (primary, unpaired). Secondary: paired by patient.
- GSE65127: score ~ IFN_c + status, lesional and non-lesional samples (primary). Secondary: paired by
  patient; correlation within healthy, non-lesional, peri-lesional and lesional groups.
- GSE53146: score ~ IFN_c + group.

**Standardisation.** Within each dataset, IFN_c and every senescence score are z-scored before
modelling, so slopes are in standard-deviation units and comparable across scores. The primary
readout is the standardised slope β for IFN_c, with 95% confidence interval and p-value.

**Classification of each validated score in each dataset:**
- **Coupled:** β > 0 and p < 0.05.
- **Not coupled:** p >= 0.05 and β < 0.5 × the SenMayo β in the same dataset.
- **Inconclusive:** anything else.

**Verdict per score:** the class that holds in at least 2 of the 3 datasets; otherwise inconclusive.

## 5. Interpretation rules (fixed)

- **R1, coupling reflects inflammatory secreted content:** SenMayo is coupled in at least 2 datasets
  and every validated non-secreted score (S1, S2, S3) is not coupled in at least 2 datasets.
- **R2, coupling is not explained by secreted factors:** at least one validated non-secreted score
  is coupled in at least 2 datasets. Follow-up with single-cell data is then needed, and the CDKN1A
  and p53-target caveat (below) is stated alongside.
- **R3, inconclusive:** anything else. Reported as inconclusive.

## 6. Sensitivity analyses and controls (reported, no verdict attached)

- Add the E2F proliferation score as a covariate to each Step B model, to check whether any coupling
  is carried by proliferation changes in inflamed or thickened skin.
- Lesional versus non-lesional score differences (limma, paired where possible).
- Report which genes were removed as shared with the IFN-γ set, and array coverage per score.

## 7. Multiplicity, power and caveats

- Step B has 3 scores × 3 datasets = 9 primary tests. Unadjusted and Benjamini-Hochberg-adjusted
  p-values are both reported. Verdicts rely on agreement across datasets (2 of 3); any claim of
  coupling (R2) is also stated with the adjusted p-values.
- Cohorts are small (5-15 patients). "Not coupled" is weak evidence of absence, so confidence
  interval widths are reported next to every classification.
- CDKN1A and GADD45A can be induced by interferon and by ordinary p53 activity, and the TANG set is
  p53-driven, so a coupling in S1 or S2 does not by itself show senescence.
- Fibroblast-derived validation may not transfer fully to epidermal or whole-skin samples.
- Tissue is bulk, so cell-type origin of any signal cannot be resolved here.

## 8. What will not be done

- No change to scripts 01-03 or their results.
- No change to the scores, rules or thresholds above after seeing vitiligo results.
- All runs are saved with their console output and session information in `results/`.

## 9. Amendment log

- 2026-10-27, before any analysis in this plan was run: GSE49863 and GSE6762 were checked as
  validation candidates and excluded. Both are two-colour ratio arrays (GSE49863: one series-matrix
  sample labelled as a ratio; GSE6762: four samples each labelled senescent versus proliferating),
  with no expression values and no separate proliferating group. Section 3 was changed to allow
  RNA-seq datasets and to list GSE130100 and GSE42509 as candidates.
- 2026-10-27, before any analysis in this plan was run: GSE42509 was checked and excluded. Its
  RNA-seq part has only one sample labelled as the fully senescent group (S.rna) and two
  pre-senescent samples, below the minimum of 3 samples per group. Its other 8 samples are
  ribosome profiling, not transcriptome data. GSE130100 passes the check on sample numbers
  (3 empty-vector controls, 3 HRAS samples); whether its supplementary table can be scored is
  still to be confirmed.
- 2026-10-28, before any analysis in this plan was run: the validation datasets are fixed as
  GSE130100 (RAS-induced senescence in BJ-TERT cells; 3 empty-vector controls versus 3 HRAS samples;
  group labels checked with HRAS expression) and GSE63577 (RNA-seq, five fibroblast lines, 3
  replicates per group; the study describes the comparison as senescent versus young proliferating
  fibroblasts). That the late-passage cells are senescent is taken from the study description and
  was not independently verified.
  GSE63577 counts as ONE validation dataset. Its primary test is score ~ group + cell line over the
  30 young versus late-passage samples (MRC-5 PD32 vs PD72, BJ young vs old, IMR-90, WI-38, HFF
  PD16 vs PD74). Per-cell-line differences and the intermediate-passage samples (HFF PD26/46/64,
  MRC-5 PD42/52/62) are reported descriptively, with no verdict. GSE130099 is not used.
  Because only two validation datasets are used, "at least 2" means a score must pass in both.
  RNA-seq processing: library-size CPM with no further normalisation, genes kept if CPM > 1 in at
  least 3 samples, log2(CPM + 1); counts of rows sharing a gene symbol are summed. Scores that are
  not validated are dropped from Step B.
  Observations made while checking the data, before any scoring: in raw counts LMNB1 was lower in
  senescent samples in both datasets, CDKN1A was not clearly higher in GSE130100, and CDKN2A changed
  little in the BJ cells of GSE63577. These do not change any rule.
