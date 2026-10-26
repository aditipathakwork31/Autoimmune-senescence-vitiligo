# Autoimmune senescence: starting with vitiligo

An exploratory, hypothesis-generating analysis of public gene-expression data.
Work in progress.

## Question

Do senescence-associated gene signatures and IFN-γ / JAK-STAT signatures rise together in
vitiligo skin, and are they separable?

## Data

All data are public (NCBI GEO). No patient-level or personal data are included in this repository.

| Dataset | Tissue | Platform | Samples | Design |
|---|---|---|---|---|
| GSE75819 | epidermis | Illumina | 15 patients, 30 arrays | paired lesional / non-lesional |
| GSE65127 | whole skin | Affymetrix GPL570 | 40 arrays | healthy (10), lesional, non-lesional, peri-lesional (10 each); lesional groups paired by patient |
| GSE53146 | FFPE skin biopsy | Illumina GPL14951 | 10 arrays | vitiligo (5) vs control (5), unpaired |

## Methods (summary)

1. Download each dataset with GEOquery; labels taken from GEO sample annotation.
2. Normalization: GSE75819 quantile-normalized and log2-transformed (the two groups of arrays
   differed in intensity; processing batch was used as a covariate); GSE65127 used as
   provided (RMA, log2); GSE53146 log2-transformed (already quantile-normalized).
3. Probes collapsed to genes (mean of probes per gene symbol).
4. Pathway scores per sample with GSVA: MSigDB Hallmark IFN-γ response, IFN-α response and
   IL6-JAK-STAT3; senescence signatures SenMayo (SAUL_SEN_MAYO), Fridman senescence (up) and
   Reactome SASP.
5. Genes shared between the IFN-γ set and any senescence set were removed before testing
   whether the scores are associated ("cleaned" sets).
6. Lesional vs non-lesional comparison with limma (paired by patient where possible);
   association between IFN-γ and senescence scores with linear models (adjusted for status,
   batch, and patient where applicable).
7. Sanity checks: CXCL9, CXCL10, STAT1 (interferon axis) and TYR, DCT, MLANA (melanocytes).

## Results so far

- Melanocyte genes (TYR, DCT, MLANA) are lower in lesional skin in all three datasets.
- IFN-γ and IFN-α response scores are higher in lesions in GSE75819 and GSE53146, but not in
  GSE65127 (whole skin).
- A SenMayo score is positively associated with the IFN-γ score in all three datasets. In
  GSE65127 the direction is also positive in healthy and non-lesional skin (n = 10 per group,
  not individually significant), so the association does not appear specific to vitiligo.
- Within-patient (paired) estimates of the SenMayo ~ IFN-γ association are mixed: not
  significant in GSE75819 (about 0.23, p about 0.30) and significant in GSE65127 (about 0.36,
  p about 0.001). Evidence for an association within the same person is therefore inconclusive.
- Senescence scores in lesions relative to normal skin are inconsistent across datasets.
- Fridman and Reactome SASP signatures do not replicate. SASP is poorly covered on the
  GSE53146 array (43 of 107 genes).

Score tables and figures are in the `results/` folders. Test statistics (estimates, standard
errors, p-values) are printed to the console when each script is run.

## Limitations

Small cohorts (5-15 patients), bulk tissue (cell types are mixed), different platforms and
tissue types, no disease-activity or treatment information, correlation only. Senescence
signatures such as SenMayo contain many inflammatory secreted factors, so they can track
interferon activity without indicating senescent cells.

## Reproduce

Requires R and the packages GEOquery, limma, GSVA, msigdbr, ggplot2. Run the scripts in
`scripts/` in order. Each script downloads its dataset from GEO and writes results to a
folder named `results_<accession>` under `~/autoimmune-senescence/` (in this repository these
were moved into `results/` and renamed to the accession number).

## Repository layout

- `scripts/` analysis scripts (one per dataset)
- `results/` figures, score tables and session information for each dataset

## Next

Single-cell skin data to ask which cell types carry a senescence signal; a senescence
signature that does not rely on secreted inflammatory factors.
