# ============================================================
# Script 04: does the IFN-gamma coupling of senescence scores survive when the
# score is built from non-secreted markers?   Follows docs/analysis_plan_04.md
#
# Step A  validate the scores in senescent fibroblasts
#         (GSE130100 RAS-induced; GSE63577 replicative, 5 cell lines)
# Step B  test the coupling to IFN-gamma in GSE75819, GSE65127, GSE53146
#         (only if SenMayo is validated in Step A, as the plan requires)
#
# Gene lists, thresholds and decision rules are fixed in the plan and in the
# CONFIG section below. Do not edit them after seeing vitiligo results.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(GSVA)
  library(msigdbr);  library(ggplot2); library(readxl)
})
options(timeout = 600)

## ---------------- CONFIG ----------------
REPO   <- path.expand("~/autoimmune-repo")      # local clone of the GitHub repository
OUTDIR <- file.path(REPO, "results", "04_nonsecreted_senescence")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (sink.number() > 0) sink()                  # close a log left open by an earlier failed run
sink(file.path(OUTDIR, "console_output.txt"), split = TRUE)

# Frozen by the plan (Section 2). Changes allowed only after Step A, logged in git.
PANEL_UP   <- c("CDKN1A", "CDKN2A", "CDKN2B", "GADD45A")
PANEL_DOWN <- c("LMNB1")
DROP_GENES <- character(0)                     # genes dropped after Step A, if any (log in git)

## ==== FUNCTIONS START ====

zs <- function(v) as.numeric(scale(v))

# counts (genes x samples) -> filter low genes -> log2 CPM  (plan, Section 3)
log_cpm <- function(counts) {
  cpm  <- t(t(counts) / colSums(counts)) * 1e6
  keep <- rowSums(cpm > 1) >= 3
  log2(cpm[keep, , drop = FALSE] + 1)
}

# sum counts of rows that share a gene symbol
collapse_counts <- function(counts, symbols) {
  ok <- !is.na(symbols) & symbols != ""
  rowsum(counts[ok, , drop = FALSE], group = symbols[ok])
}

# Arrest panel: mean of z-scored genes, down-genes sign-flipped (plan, S1)
panel_score <- function(x, up, down, min_genes = 3) {
  up <- intersect(up, rownames(x)); down <- intersect(down, rownames(x))
  if (length(up) + length(down) < min_genes) return(rep(NA_real_, ncol(x)))
  z <- function(g) t(apply(x[g, , drop = FALSE], 1, zs))
  parts <- NULL
  if (length(up))   parts <- rbind(parts,  z(up))
  if (length(down)) parts <- rbind(parts, -z(down))
  colMeans(parts)
}

# all scores for one expression matrix (log scale, genes x samples)
score_all <- function(x) {
  sets <- gs[c("S0", "S2UP", "S2DN", "S3", "E2F", "IFNC")]
  n_present <- sapply(sets, function(g) sum(g %in% rownames(x)))
  ok  <- n_present >= 10                       # plan: a GSVA set needs >= 10 present genes
  out <- data.frame(row.names = colnames(x))
  if (any(ok)) {
    gv <- gsva(gsvaParam(x, sets[ok], kcdf = "Gaussian", minSize = 10))
    for (nm in rownames(gv)) out[[nm]] <- as.numeric(gv[nm, colnames(x)])
  }
  if (all(c("S2UP", "S2DN") %in% names(out))) out$S2 <- out$S2UP - out$S2DN
  out$S1 <- panel_score(x, panel_up, panel_dn)
  attr(out, "coverage") <- data.frame(
    set = names(sets), n_in_set = sapply(sets, function(g) length(unique(g))),
    n_present = n_present, scored = ok, row.names = NULL)
  out
}

# Step A test: senescent vs proliferating, optional adjustment factor
val_test <- function(sc, group, adj = NULL) {
  do.call(rbind, lapply(c("S0", "S1", "S2", "S3"), function(s) {
    if (!(s %in% names(sc)) || all(is.na(sc[[s]]))) {
      return(data.frame(score = s, estimate = NA_real_, p = NA_real_, passes = FALSE))
    }
    d <- data.frame(y = sc[[s]], group = group)
    if (!is.null(adj)) d$adj <- adj
    m  <- lm(if (is.null(adj)) y ~ group else y ~ group + adj, data = d)
    cf <- summary(m)$coefficients["groupsenescent", ]
    data.frame(score = s, estimate = cf[[1]], p = cf[[4]],
               passes = (cf[[1]] > 0 & cf[[4]] < 0.05))
  }))
}

# per-gene change (senescent minus proliferating, log2 CPM) for the panel genes
gene_changes <- function(x, group, genes) {
  genes <- genes[genes %in% rownames(x)]
  data.frame(gene = genes, diff = sapply(genes, function(g)
    mean(x[g, group == "senescent"]) - mean(x[g, group == "proliferating"])), row.names = NULL)
}

# Step B: standardised IFN_c slope for each score, one dataset
coupling_one <- function(sc, meta, dsname, scores = c("S0", "S1", "S2", "S3")) {
  d <- cbind(sc, meta)
  d$zIFN <- zs(d$IFNC); d$zE2F <- zs(d$E2F)
  base <- switch(dsname,
    GSE75819 = "zS ~ zIFN + batch + status",
    GSE65127 = "zS ~ zIFN + status",
    GSE53146 = "zS ~ zIFN + group")
  has_patient <- "patient" %in% colnames(meta)
  rows <- lapply(scores, function(s) {
    if (!(s %in% names(d)) || all(is.na(d[[s]]))) return(NULL)
    d$zS <- zs(d[[s]])
    m   <- lm(as.formula(base), data = d)
    cf  <- summary(m)$coefficients["zIFN", ]
    ci  <- confint(m)["zIFN", ]
    me  <- lm(as.formula(paste(base, "+ zE2F")), data = d)
    cfe <- summary(me)$coefficients["zIFN", ]
    bp <- NA_real_; pp <- NA_real_
    if (has_patient) {
      mp <- lm(zS ~ zIFN + status + patient, data = d)
      cp <- summary(mp)$coefficients["zIFN", ]
      bp <- cp[[1]]; pp <- cp[[4]]
    }
    data.frame(dataset = dsname, score = s, beta = cf[[1]], ci_low = ci[[1]], ci_high = ci[[2]],
               p = cf[[4]], beta_E2Fadj = cfe[[1]], p_E2Fadj = cfe[[4]],
               beta_paired = bp, p_paired = pp)
  })
  do.call(rbind, rows)
}

# classification rule from the plan (Section 4)
classify <- function(beta, p, beta_ref) {
  ifelse(beta > 0 & p < 0.05, "coupled",
         ifelse(p >= 0.05 & beta < 0.5 * beta_ref, "not_coupled", "inconclusive"))
}

# verdict: the class that holds in at least 2 of the datasets
verdict_of <- function(cls) {
  tb <- table(cls); w <- names(tb)[tb >= 2]
  if (length(w)) w[1] else "inconclusive"
}

# score differences between disease and normal skin (secondary, no verdict)
group_diff <- function(sc, meta, dsname, scores = c("S0", "S1", "S2", "S3", "IFNC", "E2F")) {
  d <- cbind(sc, meta)
  do.call(rbind, lapply(scores, function(s) {
    if (!(s %in% names(d)) || all(is.na(d[[s]]))) return(NULL)
    m <- if (dsname == "GSE53146") lm(d[[s]] ~ group, data = d) else lm(d[[s]] ~ status + patient, data = d)
    cf <- summary(m)$coefficients
    r  <- grep("^(statuslesional|groupvitiligo)$", rownames(cf))
    data.frame(dataset = dsname, score = s, disease_minus_normal = cf[r, 1], p = cf[r, 4])
  }))
}

## ==== END FUNCTIONS ====

## ---------------- GENE SETS (frozen) ----------------
h  <- msigdbr(species = "Homo sapiens", collection = "H")
c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
allsets <- rbind(h[, c("gs_name", "gene_symbol")], c2[, c("gs_name", "gene_symbol")])
S <- split(allsets$gene_symbol, allsets$gs_name)

ifn <- S$HALLMARK_INTERFERON_GAMMA_RESPONSE
sen_old <- list(S$SAUL_SEN_MAYO, S$FRIDMAN_SENESCENCE_UP,
                S$REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP)
rm_ifn <- function(g) setdiff(g, ifn)          # genes shared with IFN-gamma are removed

gs <- list(
  S0   = rm_ifn(S$SAUL_SEN_MAYO),
  S2UP = rm_ifn(S$TANG_SENESCENCE_TP53_TARGETS_UP),
  S2DN = rm_ifn(S$TANG_SENESCENCE_TP53_TARGETS_DN),
  S3   = rm_ifn(S$CHICAS_RB1_TARGETS_SENESCENT),
  E2F  = S$HALLMARK_E2F_TARGETS,                              # proliferation control
  IFNC = setdiff(ifn, unique(unlist(sen_old))))               # same IFN_c as scripts 01-03
panel_up <- setdiff(setdiff(PANEL_UP, DROP_GENES), ifn)
panel_dn <- setdiff(setdiff(PANEL_DOWN, DROP_GENES), ifn)

cat("Genes removed from each senescence set because they are in the IFN-gamma Hallmark set:\n")
removed <- list(
  SenMayo = intersect(S$SAUL_SEN_MAYO, ifn),
  TANG_UP = intersect(S$TANG_SENESCENCE_TP53_TARGETS_UP, ifn),
  TANG_DN = intersect(S$TANG_SENESCENCE_TP53_TARGETS_DN, ifn),
  CHICAS  = intersect(S$CHICAS_RB1_TARGETS_SENESCENT, ifn),
  panel   = intersect(c(PANEL_UP, PANEL_DOWN), ifn))
print(lapply(removed, function(g) if (length(g)) paste(g, collapse = " ") else "none"))
cat("Panel genes used: up =", paste(panel_up, collapse = ", "),
    "| down =", paste(panel_dn, collapse = ", "), "\n")
cat("Set sizes:\n"); print(sapply(gs, length))

## ============================================================
## STEP A: validation in senescent fibroblasts
## ============================================================
cat("\n================ STEP A: validation in fibroblasts ================\n")

val_130100 <- function() {
  f   <- getGEOSuppFiles("GSE130100", makeDirectory = FALSE, baseDir = tempdir())
  tmp <- tempfile(fileext = ".csv")
  # keep only the first 8 comma-separated fields (gene id, symbol, 6 count columns);
  # later text fields contain commas and confuse read.csv
  system(paste("gunzip -c", shQuote(rownames(f)[1]), "| cut -d, -f1-8 >", shQuote(tmp)))
  d   <- read.csv(tmp, check.names = FALSE)
  cnt <- as.matrix(d[, 3:8]); storage.mode(cnt) <- "numeric"
  stopifnot(ncol(cnt) == 6, all(grepl("^(pBABE|HRAS)", colnames(cnt))))
  cnt <- collapse_counts(cnt, d$hgnc_symbol)
  grp <- factor(ifelse(grepl("^HRAS", colnames(cnt)), "senescent", "proliferating"),
                levels = c("proliferating", "senescent"))
  list(x = log_cpm(cnt), group = grp)
}

val_63577 <- function() {
  f  <- getGEOSuppFiles("GSE63577", makeDirectory = FALSE, baseDir = tempdir())
  gz <- rownames(f)
  for (z in gz) system2("gunzip", c("-k", "-f", shQuote(z)))
  xf <- sub("\\.gz$", "", gz)
  c1 <- as.data.frame(read_excel(grep("jenage",   xf, value = TRUE), sheet = "counts"))
  c2 <- as.data.frame(read_excel(grep("RPKMs_II", xf, value = TRUE), sheet = "counts"))
  s1 <- names(c1)[5:ncol(c1)]; s2 <- names(c2)[5:ncol(c2)]
  m  <- merge(c1[, c("ensembl_gene_id", "external_gene_id", s1)],
              c2[, c("ensembl_gene_id", s2)], by = "ensembl_gene_id")
  cnt <- as.matrix(m[, c(s1, s2)]); storage.mode(cnt) <- "numeric"
  stopifnot(ncol(cnt) == 48)
  cnt <- collapse_counts(cnt, m$external_gene_id)
  n <- colnames(cnt)
  line <- ifelse(grepl("^BJ", n), "BJ", ifelse(grepl("^IMR90", n), "IMR90",
          ifelse(grepl("^WI_38", n), "WI38", ifelse(grepl("^HFF", n), "HFF", "MRC5"))))
  young <- grepl("_Y[0-9]$", n) | grepl("^HFF_PD16_", n) | grepl("^MRC_5_PD32_", n)
  old   <- grepl("_OLD_[0-9]$", n) | grepl("_O[0-9]$", n) | grepl("^HFF_PD74_", n) | grepl("^MRC_5_PD72_", n)
  grp <- ifelse(young, "proliferating", ifelse(old, "senescent", NA))
  pd  <- suppressWarnings(as.numeric(sub(".*PD([0-9]+)_.*", "\\1", ifelse(grepl("PD[0-9]+_", n), n, NA))))
  list(x = log_cpm(cnt), group = factor(grp, levels = c("proliferating", "senescent")),
       line = line, pd = pd)
}

v1 <- val_130100(); sc1 <- score_all(v1$x)
cat("\nGSE130100: genes after filtering:", nrow(v1$x), "| samples:", ncol(v1$x), "\n")
cat("Group sizes:\n"); print(table(v1$group))
r1 <- val_test(sc1, v1$group)

v2 <- val_63577(); sc2 <- score_all(v2$x)
cat("\nGSE63577: genes after filtering:", nrow(v2$x), "| samples:", ncol(v2$x), "\n")
use <- !is.na(v2$group)
cat("Primary comparison samples (young vs late passage), by cell line:\n")
print(table(v2$line[use], v2$group[use]))
r2 <- val_test(sc2[use, ], droplevels(v2$group[use]), adj = factor(v2$line[use]))

vres <- rbind(cbind(dataset = "GSE130100", r1), cbind(dataset = "GSE63577", r2))
cat("\nValidation, senescent minus proliferating (GSE63577 adjusted for cell line):\n")
print(vres, digits = 3, row.names = FALSE)
validated <- sapply(c("S0", "S1", "S2", "S3"),
                    function(s) sum(vres$passes[vres$score == s], na.rm = TRUE) >= 2)
cat("\nValidated (positive and p < 0.05 in both datasets):\n"); print(validated)
write.csv(vres, file.path(OUTDIR, "validation_results.csv"), row.names = FALSE)
write.csv(data.frame(score = names(validated), validated = validated),
          file.path(OUTDIR, "validation_verdicts.csv"), row.names = FALSE)

# Per-gene direction of the panel genes (report only; changes need a logged decision)
gc1 <- gene_changes(v1$x, v1$group, c(panel_up, panel_dn))
gc2 <- gene_changes(v2$x[, use], droplevels(v2$group[use]), c(panel_up, panel_dn))
gchg <- merge(gc1, gc2, by = "gene", all = TRUE, suffixes = c("_GSE130100", "_GSE63577"))
cat("\nPanel genes, senescent minus proliferating (log2 CPM):\n"); print(gchg, digits = 3)
write.csv(gchg, file.path(OUTDIR, "panel_gene_changes.csv"), row.names = FALSE)
opp <- gchg$gene[sapply(seq_len(nrow(gchg)), function(i) {
  sg <- sign(c(gchg[i, 2], gchg[i, 3])); expected <- if (gchg$gene[i] %in% panel_dn) -1 else 1
  all(!is.na(sg)) && all(sg == -expected) })]
cat("Panel genes moving the wrong way in both datasets:",
    if (length(opp)) paste(opp, collapse = ", ") else "none", "\n")

# Secondary: dose-response with passage number (HFF and MRC-5, 5 passages each)
cat("\nSecondary: Spearman correlation of score with population doublings (no verdict):\n")
dose <- do.call(rbind, lapply(c("HFF", "MRC5"), function(l) {
  i <- v2$line == l & !is.na(v2$pd)
  do.call(rbind, lapply(c("S0", "S1", "S2", "S3"), function(s) {
    if (!(s %in% names(sc2)) || all(is.na(sc2[[s]]))) return(NULL)
    data.frame(line = l, score = s, n = sum(i), rho = suppressWarnings(cor(sc2[[s]][i], v2$pd[i], method = "spearman")))
  }))
}))
print(dose, digits = 3, row.names = FALSE)
write.csv(dose, file.path(OUTDIR, "validation_dose_response.csv"), row.names = FALSE)

# Per-cell-line descriptive differences in GSE63577
perline <- do.call(rbind, lapply(unique(v2$line), function(l) {
  i <- use & v2$line == l
  do.call(rbind, lapply(c("S0", "S1", "S2", "S3"), function(s) {
    if (!(s %in% names(sc2)) || all(is.na(sc2[[s]]))) return(NULL)
    data.frame(line = l, score = s,
               diff = mean(sc2[[s]][i & v2$group == "senescent"], na.rm = TRUE) -
                      mean(sc2[[s]][i & v2$group == "proliferating"], na.rm = TRUE))
  }))
}))
cat("\nGSE63577 per cell line, senescent minus proliferating (descriptive):\n")
print(reshape(perline, idvar = "line", timevar = "score", direction = "wide"), digits = 3, row.names = FALSE)

# Figure: validation
lg <- function(sc, group, ds) do.call(rbind, lapply(c("S0", "S1", "S2", "S3"), function(s)
  if (s %in% names(sc)) data.frame(dataset = ds, score = s, group = group, value = sc[[s]])))
vl <- rbind(lg(sc1, v1$group, "GSE130100"), lg(sc2[use, ], droplevels(v2$group[use]), "GSE63577"))
pv <- ggplot(vl, aes(group, value, color = group)) + geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 1.5) + facet_grid(score ~ dataset, scales = "free") +
  labs(title = "Step A: scores in proliferating vs senescent fibroblasts", x = NULL, y = "score") +
  theme(legend.position = "none")
ggsave(file.path(OUTDIR, "validation_scores.png"), pv, width = 6.5, height = 8, dpi = 200)

STEP_B <- isTRUE(validated[["S0"]])
if (!STEP_B) {
  cat("\n*** STOP: SenMayo (S0) was not validated in both fibroblast datasets. ",
      "The plan says to stop and revisit it before Step B. Step B was not run. ***\n", sep = "")
}

## ============================================================
## STEP B: coupling to IFN-gamma in the vitiligo datasets
## ============================================================
if (STEP_B) {
cat("\n================ STEP B: coupling in vitiligo datasets ================\n")

load_75819 <- function() {   # same preprocessing as script 01
  e  <- getGEO("GSE75819", GSEMatrix = TRUE)[[1]]; ph <- pData(e)
  status  <- ifelse(ph$`tissue:ch1` == "Lesional Vitiligo Skin", "lesional", "nonlesional")
  patient <- sub("[NV]$", "", ph$title)
  x <- exprs(e); batch <- ifelse(apply(x, 2, median) < 20, "A", "B")
  xl <- log2(pmax(normalizeBetweenArrays(x, method = "quantile"), 1))
  sym <- as.character(fData(e)$ILMN_Gene); keep <- !is.na(sym) & sym != ""
  xg <- avereps(xl[keep, ], ID = sym[keep])
  list(xg = xg, meta = data.frame(status = factor(status, levels = c("nonlesional", "lesional")),
       patient = factor(patient), batch = batch, row.names = colnames(xg)))
}
load_65127 <- function() {   # same preprocessing as script 02
  e  <- getGEO("GSE65127", GSEMatrix = TRUE)[[1]]; ph <- pData(e)
  grp <- c("Healthy Volunteer" = "healthy", "Lesional" = "lesional",
           "Non Lesional" = "nonlesional", "Peri-Lesional" = "perilesional")[as.character(ph$source_name_ch1)]
  patient <- sub("^.*\\[([0-9]+)_S[0-9]+_[A-Z]+\\]$", "\\1", ph$title)
  x <- exprs(e); if (max(x) > 100) x <- log2(pmax(x, 1))
  sym <- sub(" ///.*$", "", as.character(fData(e)[["Gene Symbol"]]))
  keep <- !is.na(sym) & sym != "" & sym != "---"
  xg <- avereps(x[keep, ], ID = sym[keep])
  list(xg = xg, meta = data.frame(group = grp, status = factor(grp, levels = c("nonlesional", "lesional")),
       patient = factor(patient), row.names = colnames(xg)))
}
load_53146 <- function() {   # same preprocessing as script 03
  e  <- getGEO("GSE53146", GSEMatrix = TRUE)[[1]]; ph <- pData(e)
  grp <- factor(ifelse(grepl("vitiligo", ph$source_name_ch1, ignore.case = TRUE), "vitiligo", "control"),
                levels = c("control", "vitiligo"))
  x <- exprs(e); if (max(x) > 100) x <- log2(pmax(x, 1))
  sym <- as.character(fData(e)[["ILMN_Gene"]]); keep <- !is.na(sym) & sym != ""
  xg <- avereps(x[keep, ], ID = sym[keep])
  list(xg = xg, meta = data.frame(group = grp, row.names = colnames(xg)))
}

dsl <- list(GSE75819 = load_75819(), GSE65127 = load_65127(), GSE53146 = load_53146())
res <- list(); gd <- list(); cov_all <- list()
non_val     <- intersect(c("S1", "S2", "S3"), names(validated)[validated])
keep_scores <- c("S0", non_val)
dropped     <- setdiff(c("S1", "S2", "S3"), non_val)
cat("Scores carried into Step B:", paste(keep_scores, collapse = ", "), "\n")
cat("Dropped because they were not validated in Step A:",
    if (length(dropped)) paste(dropped, collapse = ", ") else "none", "\n")
for (nm in names(dsl)) {
  cat("\n--", nm, "--\n")
  sc   <- score_all(dsl[[nm]]$xg)
  meta <- dsl[[nm]]$meta
  cov_all[[nm]] <- cbind(dataset = nm, attr(sc, "coverage"))
  # reproduction check against the earlier results (IFN_c and SenMayo scores)
  old <- tryCatch(read.csv(file.path(REPO, "results", nm, "scores_clean.csv"), row.names = 1),
                  error = function(e) NULL)
  if (!is.null(old) && all(rownames(sc) %in% rownames(old))) {
    cat("Check vs earlier scores: cor(IFN_c) =", round(cor(sc$IFNC, old[rownames(sc), "IFN_c"]), 4),
        "| cor(SenMayo) =", round(cor(sc$S0, old[rownames(sc), "SENMAYO"]), 4), "\n")
  }
  if (nm == "GSE65127") {           # primary models use lesional and non-lesional samples only
    ln <- meta$group %in% c("lesional", "nonlesional")
    sc <- sc[ln, ]; meta <- meta[ln, ]; meta$patient <- droplevels(meta$patient)
  }
  res[[nm]] <- coupling_one(sc, meta, nm, scores = keep_scores)
  gd[[nm]]  <- group_diff(sc, meta, nm)
}

cp <- do.call(rbind, res); rownames(cp) <- NULL
ref <- cp[cp$score == "S0", c("dataset", "beta")]; names(ref)[2] <- "beta_S0"
cp <- merge(cp, ref, by = "dataset")
cp$class <- classify(cp$beta, cp$p, cp$beta_S0)
prim <- cp$score %in% non_val
cp$p_BH <- NA_real_; cp$p_BH[prim] <- p.adjust(cp$p[prim], method = "BH")
cp <- cp[order(cp$score, cp$dataset), ]
cat("\nPrimary results (standardised slope of IFN_c; CI = 95%):\n")
print(cp[, c("dataset", "score", "beta", "ci_low", "ci_high", "p", "p_BH", "class")], digits = 3, row.names = FALSE)
cat("\nSensitivity: adjusting for the E2F proliferation score, and paired by patient:\n")
print(cp[, c("dataset", "score", "beta_E2Fadj", "p_E2Fadj", "beta_paired", "p_paired")], digits = 3, row.names = FALSE)
write.csv(cp, file.path(OUTDIR, "coupling_results.csv"), row.names = FALSE)

verdicts <- sapply(keep_scores, function(s) verdict_of(cp$class[cp$score == s]))
cat("\nVerdict per score (class that holds in >= 2 of 3 datasets):\n"); print(verdicts)

rule <- if (length(non_val) && any(verdicts[non_val] == "coupled")) "R2"
        else if (verdicts[["S0"]] == "coupled" && length(non_val) && all(verdicts[non_val] == "not_coupled")) "R1"
        else "R3"
cat("\nInterpretation rule triggered:", rule, "\n")
cat(switch(rule,
  R1 = "R1: coupling reflects inflammatory secreted content (SenMayo coupled, validated non-secreted scores not coupled).\n",
  R2 = "R2: coupling is not explained by secreted factors (a validated non-secreted score is coupled in >= 2 datasets). Single-cell follow-up needed; see plan caveats.\n",
  R3 = "R3: inconclusive.\n"))

gdf <- do.call(rbind, gd); rownames(gdf) <- NULL
cat("\nSecondary: score differences, disease minus normal skin (no verdict):\n"); print(gdf, digits = 3, row.names = FALSE)
write.csv(gdf, file.path(OUTDIR, "disease_vs_normal_scores.csv"), row.names = FALSE)
covdf <- do.call(rbind, cov_all); rownames(covdf) <- NULL
cat("\nGene-set coverage per dataset:\n"); print(covdf, row.names = FALSE)
write.csv(covdf, file.path(OUTDIR, "set_coverage.csv"), row.names = FALSE)
write.csv(data.frame(score = names(verdicts), verdict = verdicts, rule = rule),
          file.path(OUTDIR, "verdicts.csv"), row.names = FALSE)

pf <- ggplot(cp, aes(beta, score, color = dataset)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(width = 0.5)) +
  labs(title = "Step B: standardised slope of senescence score on IFN_c", x = "standardised slope (95% CI)", y = NULL)
ggsave(file.path(OUTDIR, "coupling_forest.png"), pf, width = 7, height = 4.5, dpi = 200)
}  # end Step B

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("\nDone. Results saved in", OUTDIR, "\n")
sink()
