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
REPO     <- path.expand("~/autoimmune-repo")   # local clone of the GitHub repository
OUTDIR   <- file.path(REPO, "results", ACC)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (sink.number() > 0) sink()                  # close a log left open by an earlier failed run
sink(file.path(OUTDIR, "console_output.txt"), split = TRUE)   # copy of everything printed

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
ifn_chk <- chk(c("CXCL9", "CXCL10", "STAT1", "IFNG", "GZMB"))
mel_chk <- chk(c("MITF", "TYR", "DCT", "PMEL", "MLANA"))
cat("\nInterferon genes, lesional minus non-lesional (log2):\n"); print(ifn_chk)
cat("\nMelanocyte genes, lesional minus non-lesional (log2):\n"); print(mel_chk)
write.csv(data.frame(gene = c(names(ifn_chk), names(mel_chk)),
                     group = c(rep("interferon", length(ifn_chk)), rep("melanocyte", length(mel_chk))),
                     log2_diff = c(ifn_chk, mel_chk)),
          file.path(OUTDIR, "sanity_genes.csv"), row.names = FALSE)

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
tab <- topTable(fit, coef = "stlesional", number = Inf)
cat("\nPaired lesional vs non-lesional, pathway scores:\n")
print(round(tab, 3))
write.csv(tab, file.path(OUTDIR, "paired_lesional_vs_nonlesional.csv"))

## ---------------- 8. IFN vs SENESCENCE ----------------
d <- data.frame(t(sc_clean[, ln]), status = st, patient = pt)
assoc <- do.call(rbind, lapply(c("SENMAYO", "FRIDMAN", "SASP"), function(s) {
  m1 <- summary(lm(d[[s]] ~ IFN_c + status, data = d))$coefficients["IFN_c", ]
  m2 <- summary(lm(d[[s]] ~ IFN_c + status + patient, data = d))$coefficients["IFN_c", ]
  data.frame(senescence_set = s,
             model = c("unpaired_adj_status", "paired_by_patient"),
             rbind(m1, m2), row.names = NULL)
}))
names(assoc)[3:6] <- c("estimate", "std_error", "t", "p")
cat("\nIFN_c -> senescence score (lesional and non-lesional samples only):\n")
print(assoc, digits = 3)
write.csv(assoc, file.path(OUTDIR, "ifn_vs_senescence_models.csv"), row.names = FALSE)

# Control: the same association within each group, including healthy volunteers
dd <- data.frame(t(sc_clean), group = group)
ctrl <- do.call(rbind, lapply(c("healthy", "nonlesional", "perilesional", "lesional"), function(g) {
  s  <- dd[dd$group == g, ]
  ct <- suppressWarnings(cor.test(s$IFN_c, s$SENMAYO, method = "spearman"))
  data.frame(group = g, n = nrow(s), rho = unname(ct$estimate), p = ct$p.value)
}))
cat("\nSpearman correlation IFN_c vs SenMayo within each group (n = 10 each):\n")
print(ctrl, digits = 3)
write.csv(ctrl, file.path(OUTDIR, "spearman_ifn_vs_senmayo_by_group.csv"), row.names = FALSE)

pooled <- summary(lm(SENMAYO ~ IFN_c + group, data = dd))$coefficients["IFN_c", ]
cat("\nPooled model SENMAYO ~ IFN_c + group (all 40 samples; samples from the same patient are not independent):\n")
print(round(pooled, 4))

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

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("\nDone. Results saved in", OUTDIR, "\n")
sink()
