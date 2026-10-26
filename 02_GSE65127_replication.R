# ============================================================
# Replication of the senescence vs IFN-gamma analysis
# Dataset: GSE65127 (vitiligo whole skin, Affymetrix GPL570, RMA)
# Groups: healthy (10), lesional (10), non-lesional (10), peri-lesional (10)
# Lesional / non-lesional / peri-lesional samples are paired by patient.
# To reuse for another dataset, edit only the CONFIG and LABELS sections.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(GSVA)
  library(msigdbr);  library(ggplot2)
})
options(timeout = 600)

## ---------------- CONFIG ----------------
ACC      <- "GSE65127"
GENE_COL <- "Gene Symbol"          # column in fData() holding gene symbols
OUTDIR   <- "~/autoimmune-senescence/results_GSE65127"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## ---------------- 1. DOWNLOAD ----------------
gset  <- getGEO(ACC, GSEMatrix = TRUE)
eset  <- gset[[1]]
pheno <- pData(eset)

## ---------------- 2. LABELS (dataset-specific) ----------------
group <- c("Healthy Volunteer" = "healthy",
           "Lesional"          = "lesional",
           "Non Lesional"      = "nonlesional",
           "Peri-Lesional"     = "perilesional")[as.character(pheno$source_name_ch1)]
patient <- sub("^.*\\[([0-9]+)_S[0-9]+_[A-Z]+\\]$", "\\1", pheno$title)

print(table(group, useNA = "ifany"))
ln <- group %in% c("lesional", "nonlesional")
stopifnot(all(table(patient[ln]) == 2))   # every patient has one of each

## ---------------- 3. EXPRESSION ----------------
x <- exprs(eset)
if (max(x) > 100) x <- log2(pmax(x, 1))   # RMA data is already log2
cat("Value range:", range(x), "\n")

sym  <- as.character(fData(eset)[[GENE_COL]])
sym  <- sub(" ///.*$", "", sym)           # keep first symbol if a probe maps to several
keep <- !is.na(sym) & sym != "" & sym != "---"
xg   <- avereps(x[keep, ], ID = sym[keep])
cat("Genes x samples:", dim(xg), "\n")

## ---------------- 4. SANITY CHECKS ----------------
chk <- function(genes, a = "lesional", b = "nonlesional") {
  genes <- genes[genes %in% rownames(xg)]
  round(sapply(genes, function(g)
    mean(xg[g, group == a]) - mean(xg[g, group == b])), 2)
}
cat("\nInterferon genes, lesional minus non-lesional (log2):\n")
print(chk(c("CXCL9", "CXCL10", "STAT1", "IFNG", "GZMB")))
cat("\nMelanocyte genes, lesional minus non-lesional (log2):\n")
print(chk(c("MITF", "TYR", "DCT", "PMEL", "MLANA")))

## ---------------- 5. GENE SETS ----------------
h  <- msigdbr(species = "Homo sapiens", collection = "H")
c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
allsets <- rbind(h[, c("gs_name", "gene_symbol")], c2[, c("gs_name", "gene_symbol")])
S <- split(allsets$gene_symbol, allsets$gs_name)

ifn  <- S$HALLMARK_INTERFERON_GAMMA_RESPONSE
sen  <- list(SENMAYO = S$SAUL_SEN_MAYO,
             FRIDMAN = S$FRIDMAN_SENESCENCE_UP,
             SASP    = S$REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP)

# (a) full sets, for the lesional vs non-lesional comparison
sets_full <- list(
  IFN_GAMMA = ifn,
  IFN_ALPHA = S$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  IL6_JAK_STAT3 = S$HALLMARK_IL6_JAK_STAT3_SIGNALING,
  SENMAYO = sen$SENMAYO, FRIDMAN = sen$FRIDMAN, SASP = sen$SASP)

# (b) cleaned sets: genes shared between IFN-gamma and any senescence list removed
sets_clean <- list(
  IFN_c   = setdiff(ifn, unique(unlist(sen))),
  SENMAYO = setdiff(sen$SENMAYO, ifn),
  FRIDMAN = setdiff(sen$FRIDMAN, ifn),
  SASP    = setdiff(sen$SASP, ifn))

## ---------------- 6. SCORING ----------------
sc_full  <- gsva(gsvaParam(xg, sets_full,  kcdf = "Gaussian"))
sc_clean <- gsva(gsvaParam(xg, sets_clean, kcdf = "Gaussian"))
write.csv(t(sc_full),  file.path(OUTDIR, "scores_full.csv"))
write.csv(t(sc_clean), file.path(OUTDIR, "scores_clean.csv"))

## ---------------- 7. LESIONAL vs NON-LESIONAL (paired) ----------------
pt  <- factor(patient[ln])
st  <- factor(group[ln], levels = c("nonlesional", "lesional"))
fit <- eBayes(lmFit(sc_full[, ln], model.matrix(~ pt + st)))
cat("\nPaired lesional vs non-lesional, pathway scores:\n")
print(round(topTable(fit, coef = "stlesional", number = Inf), 3))

## ---------------- 8. IFN vs SENESCENCE ----------------
d <- data.frame(t(sc_clean[, ln]), status = st, patient = pt)
cat("\nIFN_c -> senescence score, adjusting for status (unpaired):\n")
for (s in c("SENMAYO", "FRIDMAN", "SASP")) {
  cf <- summary(lm(d[[s]] ~ IFN_c + status, data = d))$coefficients
  cat(s, ": ", paste(round(cf["IFN_c", ], 3), collapse = "  "), "\n")
}
cat("\nSame, paired by patient (patient as fixed effect):\n")
for (s in c("SENMAYO", "FRIDMAN", "SASP")) {
  cf <- summary(lm(d[[s]] ~ IFN_c + status + patient, data = d))$coefficients
  cat(s, ": ", paste(round(cf["IFN_c", ], 3), collapse = "  "), "\n")
}
cat("(columns: Estimate, Std.Error, t, p)\n")

## ---------------- 9. FIGURES ----------------
dd <- data.frame(t(sc_clean), group = group)
p1 <- ggplot(dd, aes(IFN_c, SENMAYO, color = group)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = ACC, x = "IFN-gamma score (cleaned)", y = "SenMayo score (cleaned)")
ggsave(file.path(OUTDIR, "scatter_IFN_vs_SenMayo.png"), p1, width = 6, height = 5, dpi = 200)

pca <- prcomp(t(xg[apply(xg, 1, var) > quantile(apply(xg, 1, var), 0.5), ]))
p2 <- ggplot(data.frame(pca$x[, 1:2], group = group), aes(PC1, PC2, color = group)) +
  geom_point(size = 3) + labs(title = paste(ACC, "PCA"))
ggsave(file.path(OUTDIR, "pca.png"), p2, width = 6, height = 5, dpi = 200)

cat("\nDone. Results saved in", OUTDIR, "\n")
