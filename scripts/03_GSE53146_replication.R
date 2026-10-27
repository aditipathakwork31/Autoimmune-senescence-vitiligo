# ============================================================
# Consistency check of the senescence vs IFN-gamma analysis
# Dataset: GSE53146 (FFPE skin biopsies, Illumina GPL14951)
# 5 vitiligo vs 5 control, NOT paired (different people).
# Tiny sample: use for direction only, not for significance.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(GSVA)
  library(msigdbr);  library(ggplot2)
})
options(timeout = 600)

## ---------------- CONFIG ----------------
ACC      <- "GSE53146"
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
group <- ifelse(grepl("vitiligo", pheno$source_name_ch1, ignore.case = TRUE),
                "vitiligo", "control")
group <- factor(group, levels = c("control", "vitiligo"))
print(table(group))
stopifnot(all(table(group) == 5))

## ---------------- 3. EXPRESSION ----------------
x <- exprs(eset)
if (max(x) > 100) x <- log2(pmax(x, 1))    # raw-scale values -> log2
cat("Value range after log2:", range(x), "\n")

sym  <- as.character(fData(eset)[[GENE_COL]])
keep <- !is.na(sym) & sym != ""
xg   <- avereps(x[keep, ], ID = sym[keep])
cat("Genes x samples:", dim(xg), "\n")

## ---------------- 4. SANITY CHECKS ----------------
chk <- function(genes) {
  present <- genes[genes %in% rownames(xg)]
  missing <- setdiff(genes, present)
  if (length(missing)) cat("Not on array / filtered out:", missing, "\n")
  round(sapply(present, function(g)
    mean(xg[g, group == "vitiligo"]) - mean(xg[g, group == "control"])), 2)
}
ifn_chk <- chk(c("CXCL9", "CXCL10", "STAT1", "IFNG", "GZMB"))
mel_chk <- chk(c("MITF", "TYR", "DCT", "PMEL", "MLANA"))
cat("\nInterferon genes, vitiligo minus control (log2):\n"); print(ifn_chk)
cat("\nMelanocyte genes, vitiligo minus control (log2):\n"); print(mel_chk)
write.csv(data.frame(gene = c(names(ifn_chk), names(mel_chk)),
                     group = c(rep("interferon", length(ifn_chk)), rep("melanocyte", length(mel_chk))),
                     log2_diff = c(ifn_chk, mel_chk)),
          file.path(OUTDIR, "sanity_genes.csv"), row.names = FALSE)

## ---------------- 5. GENE SETS ----------------
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

sets_clean <- list(
  IFN_c   = setdiff(ifn, unique(unlist(sen))),
  SENMAYO = setdiff(sen$SENMAYO, ifn),
  FRIDMAN = setdiff(sen$FRIDMAN, ifn),
  SASP    = setdiff(sen$SASP, ifn))

# how many genes of each set are actually on this (pre-filtered) array?
cover <- function(L) data.frame(
  set = names(L), in_set = sapply(L, function(g) length(unique(g))),
  on_array = sapply(L, function(g) sum(unique(g) %in% rownames(xg))))
cov_tab <- cover(sets_clean)
cat("\nGene-set coverage on this array:\n")
print(cov_tab, row.names = FALSE)
write.csv(cov_tab, file.path(OUTDIR, "gene_set_coverage.csv"), row.names = FALSE)

## ---------------- 6. SCORING ----------------
sc_full  <- gsva(gsvaParam(xg, sets_full,  kcdf = "Gaussian", minSize = 5))
sc_clean <- gsva(gsvaParam(xg, sets_clean, kcdf = "Gaussian", minSize = 5))
write.csv(t(sc_full),  file.path(OUTDIR, "scores_full.csv"))
write.csv(t(sc_clean), file.path(OUTDIR, "scores_clean.csv"))

## ---------------- 7. VITILIGO vs CONTROL ----------------
fit <- eBayes(lmFit(sc_full, model.matrix(~ group)))
tab <- topTable(fit, coef = "groupvitiligo", number = Inf)
cat("\nVitiligo vs control, pathway scores (unpaired, n = 5 vs 5):\n")
print(round(tab, 3))
write.csv(tab, file.path(OUTDIR, "vitiligo_vs_control.csv"))

## ---------------- 8. IFN vs SENESCENCE ----------------
d <- data.frame(t(sc_clean), group = group)
assoc <- do.call(rbind, lapply(intersect(c("SENMAYO", "FRIDMAN", "SASP"), colnames(d)), function(s) {
  cf <- summary(lm(d[[s]] ~ IFN_c + group, data = d))$coefficients["IFN_c", ]
  data.frame(senescence_set = s, estimate = cf[[1]], std_error = cf[[2]], t = cf[[3]], p = cf[[4]],
             spearman_rho_all10 = cor(d$IFN_c, d[[s]], method = "spearman"))
}))
cat("\nIFN_c -> senescence score, adjusting for group (n = 10):\n")
print(assoc, digits = 3, row.names = FALSE)
write.csv(assoc, file.path(OUTDIR, "ifn_vs_senescence_models.csv"), row.names = FALSE)

## ---------------- 9. FIGURE ----------------
p1 <- ggplot(d, aes(IFN_c, SENMAYO, color = group)) +
  geom_point(size = 3) + geom_smooth(method = "lm", se = FALSE) +
  labs(title = ACC, x = "IFN-gamma score (cleaned)", y = "SenMayo score (cleaned)")
ggsave(file.path(OUTDIR, "scatter_IFN_vs_SenMayo.png"), p1, width = 6, height = 5, dpi = 200)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("\nDone. Results saved in", OUTDIR, "\n")
sink()
