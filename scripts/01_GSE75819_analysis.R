# ============================================================
# Senescence vs IFN-gamma signatures in vitiligo epidermis
# Dataset: GSE75819 (Illumina, 15 patients, paired lesional / non-lesional)
# This script reproduces, in one run, the analysis first done step by step.
# Sample labels: tissue:ch1 in GEO annotation. N = normal skin, V = lesional.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(GSVA)
  library(msigdbr);  library(ggplot2)
})
options(timeout = 600)

## ---------------- CONFIG ----------------
ACC      <- "GSE75819"
GENE_COL <- "ILMN_Gene"
REPO     <- path.expand("~/autoimmune-repo")   # local clone of the GitHub repository
OUTDIR   <- file.path(REPO, "results", ACC)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (sink.number() > 0) sink()                  # close a log left open by an earlier failed run
sink(file.path(OUTDIR, "console_output.txt"), split = TRUE)   # copy of everything printed

## ---------------- 1. DOWNLOAD ----------------
gset  <- getGEO(ACC, GSEMatrix = TRUE)
eset  <- gset[[1]]
pheno <- pData(eset)

## ---------------- 2. LABELS ----------------
status  <- ifelse(pheno$`tissue:ch1` == "Lesional Vitiligo Skin", "lesional", "nonlesional")
patient <- sub("[NV]$", "", pheno$title)
print(table(status))
stopifnot(all(table(patient) == 2), sum(status == "lesional") == 15)

## ---------------- 3. NORMALIZATION ----------------
x <- exprs(eset)                                  # raw-scale values (BeadStudio)
cat("Raw range:", range(x), "\n")

# Two groups of arrays differ ~8-fold in median intensity (processing batch).
# The batch is defined from raw medians and used as a covariate later.
batch <- ifelse(apply(x, 2, median) < 20, "A", "B")
print(table(batch))

xq <- normalizeBetweenArrays(x, method = "quantile")
xl <- log2(pmax(xq, 1))                           # floor at 1 before log2
cat("Range after quantile + log2:", range(xl), "\n")

png(file.path(OUTDIR, "boxplot_after_normalization.png"), width = 1200, height = 700)
boxplot(xl, las = 2, outline = FALSE, main = paste(ACC, "after quantile + log2"))
dev.off()

## ---------------- 4. PROBES -> GENES ----------------
sym  <- as.character(fData(eset)[[GENE_COL]])
keep <- !is.na(sym) & sym != ""
xg   <- avereps(xl[keep, ], ID = sym[keep])
cat("Genes x samples:", dim(xg), "\n")

## ---------------- 5. SANITY CHECKS ----------------
les <- status == "lesional"
chk <- function(genes) {
  genes <- genes[genes %in% rownames(xg)]
  round(sapply(genes, function(g) mean(xg[g, les]) - mean(xg[g, !les])), 2)
}
ifn_chk <- chk(c("CXCL9", "CXCL10", "STAT1", "IFNG", "GZMB"))
mel_chk <- chk(c("MITF", "TYR", "DCT", "PMEL", "MLANA"))
cat("\nInterferon genes, lesional minus non-lesional (log2):\n"); print(ifn_chk)
cat("\nMelanocyte genes, lesional minus non-lesional (log2):\n"); print(mel_chk)
write.csv(data.frame(gene = c(names(ifn_chk), names(mel_chk)),
                     group = c(rep("interferon", length(ifn_chk)), rep("melanocyte", length(mel_chk))),
                     log2_diff = c(ifn_chk, mel_chk)),
          file.path(OUTDIR, "sanity_genes.csv"), row.names = FALSE)

## ---------------- 6. GENE SETS ----------------
h  <- msigdbr(species = "Homo sapiens", collection = "H")
c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
allsets <- rbind(h[, c("gs_name", "gene_symbol")], c2[, c("gs_name", "gene_symbol")])
S <- split(allsets$gene_symbol, allsets$gs_name)

ifn <- S$HALLMARK_INTERFERON_GAMMA_RESPONSE
sen <- list(SENMAYO = S$SAUL_SEN_MAYO,
            FRIDMAN = S$FRIDMAN_SENESCENCE_UP,
            SASP    = S$REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP)

sets_full <- list(
  IFN_GAMMA = ifn,
  IFN_ALPHA = S$HALLMARK_INTERFERON_ALPHA_RESPONSE,
  IL6_JAK_STAT3 = S$HALLMARK_IL6_JAK_STAT3_SIGNALING,
  SENMAYO = sen$SENMAYO, FRIDMAN = sen$FRIDMAN, SASP = sen$SASP)

# genes shared between the IFN-gamma set and each senescence set are removed
cat("\nGenes shared with the IFN-gamma set:\n")
print(sapply(sen, function(g) length(intersect(g, ifn))))

sets_clean <- list(
  IFN_c   = setdiff(ifn, unique(unlist(sen))),
  SENMAYO = setdiff(sen$SENMAYO, ifn),
  FRIDMAN = setdiff(sen$FRIDMAN, ifn),
  SASP    = setdiff(sen$SASP, ifn))

## ---------------- 7. SCORING ----------------
sc_full  <- gsva(gsvaParam(xg, sets_full,  kcdf = "Gaussian"))
sc_clean <- gsva(gsvaParam(xg, sets_clean, kcdf = "Gaussian"))
write.csv(t(sc_full),  file.path(OUTDIR, "scores_full.csv"))
write.csv(t(sc_clean), file.path(OUTDIR, "scores_clean.csv"))

## ---------------- 8. LESIONAL vs NON-LESIONAL (paired) ----------------
pt  <- factor(patient)
st  <- factor(status, levels = c("nonlesional", "lesional"))
fit <- eBayes(lmFit(sc_full, model.matrix(~ pt + st)))
tab <- topTable(fit, coef = "stlesional", number = Inf)
cat("\nPaired lesional vs non-lesional, pathway scores:\n")
print(round(tab, 3))
write.csv(tab, file.path(OUTDIR, "paired_lesional_vs_nonlesional.csv"))

## ---------------- 9. IFN vs SENESCENCE ----------------
d <- data.frame(t(sc_clean), status = st, patient = pt, batch = batch)

assoc <- do.call(rbind, lapply(c("SENMAYO", "FRIDMAN", "SASP"), function(s) {
  m1 <- summary(lm(d[[s]] ~ IFN_c + batch + status, data = d))$coefficients["IFN_c", ]
  m2 <- summary(lm(d[[s]] ~ IFN_c + status + patient, data = d))$coefficients["IFN_c", ]
  data.frame(senescence_set = s,
             model = c("unpaired_adj_batch_status", "paired_by_patient"),
             rbind(m1, m2), row.names = NULL)
}))
names(assoc)[3:6] <- c("estimate", "std_error", "t", "p")
cat("\nIFN_c -> senescence score (unpaired adjusts for batch and status; paired uses patient):\n")
print(assoc, digits = 3)
write.csv(assoc, file.path(OUTDIR, "ifn_vs_senescence_models.csv"), row.names = FALSE)

cors <- do.call(rbind, lapply(c("lesional", "nonlesional"), function(g) {
  s  <- d[d$status == g, ]
  ct <- suppressWarnings(cor.test(s$IFN_c, s$SENMAYO, method = "spearman"))
  data.frame(group = g, rho = unname(ct$estimate), p = ct$p.value)
}))
cat("\nSpearman correlation IFN_c vs SenMayo within each skin type:\n")
print(cors, digits = 3)
write.csv(cors, file.path(OUTDIR, "spearman_ifn_vs_senmayo_by_status.csv"), row.names = FALSE)

## ---------------- 10. FIGURES ----------------
p1 <- ggplot(d, aes(IFN_c, SENMAYO, color = status)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = ACC, x = "IFN-gamma score (cleaned)", y = "SenMayo score (cleaned)")
ggsave(file.path(OUTDIR, "scatter_IFN_vs_SenMayo.png"), p1, width = 6, height = 5, dpi = 200)

p2 <- ggplot(d, aes(IFN_c, SENMAYO, color = batch, shape = status)) +
  geom_point(size = 3) +
  labs(title = paste(ACC, "- batch check"), x = "IFN-gamma score (cleaned)", y = "SenMayo score (cleaned)")
ggsave(file.path(OUTDIR, "scatter_batch_check.png"), p2, width = 6.5, height = 5, dpi = 200)

pca <- prcomp(t(xl))
p3 <- ggplot(data.frame(pca$x[, 1:2], batch = batch, status = status),
             aes(PC1, PC2, color = batch, shape = status)) +
  geom_point(size = 3) + labs(title = paste(ACC, "PCA"))
ggsave(file.path(OUTDIR, "pca.png"), p3, width = 6.5, height = 5, dpi = 200)

## ---------------- 11. RECORD THE SOFTWARE VERSIONS ----------------
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("\nDone. Results saved in", OUTDIR, "\n")
sink()
