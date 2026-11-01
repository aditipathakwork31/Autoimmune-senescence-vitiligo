# ============================================================
# Script 05: which skin cell types carry the interferon-senescence coupling?
# Follows docs/analysis_plan_05.md (Parts 1 and 2, dataset GSE298871).
# Part 3 (GSE203262) is not run: its lesional / non-lesional barcode mapping is undocumented.
#
# Input: the 10x files of GSE298871, extracted to DATA_DIR (outside the repository):
#   tar -xf GSE298871_RAW.tar -C ~/sc_check/inspect2
# All rules, markers and thresholds are fixed in the plan. Do not edit them after seeing results.
# ============================================================

suppressPackageStartupMessages({
  library(GEOquery); library(GSVA); library(msigdbr); library(ggplot2); library(Matrix)
})
options(timeout = 600)

## ---------------- CONFIG ----------------
REPO     <- path.expand("~/autoimmune-repo")
DATA_DIR <- path.expand("~/sc_check/inspect2")
OUTDIR   <- file.path(REPO, "results", "05_singlecell_vitiligo")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
if (sink.number() > 0) sink()                  # close a log left open by an earlier failed run
sink(file.path(OUTDIR, "console_output.txt"), split = TRUE)

## ==== FUNCTIONS START ====

MARKERS <- list(
  keratinocyte = c("KRT5", "KRT14", "KRT15", "KRT1", "KRT10", "KRT17", "KRTDAP", "SFN"),
  T_cell       = c("CD3D", "CD3E", "CD2", "TRAC", "TRBC2", "LCK"),
  myeloid      = c("LYZ", "CD68", "AIF1", "CD14", "FCER1G", "TYROBP"),
  fibroblast   = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA"),
  endothelial  = c("PECAM1", "VWF", "CLDN5", "EGFL7", "PLVAP", "ACKR1"),
  melanocyte   = c("PMEL", "MLANA", "TYR", "DCT", "MITF", "TYRP1"),
  mast         = c("TPSAB1", "TPSB2", "CPA3", "MS4A2"),
  mural        = c("ACTA2", "TAGLN", "RGS5", "MYL9"))
CLASSES     <- c(names(MARKERS), "unassigned")
PRIMARY     <- c("keratinocyte", "T_cell", "myeloid", "fibroblast")
MIN_CELLS   <- 50      # cells needed for a (sample, class) unit to be eligible
MIN_SAMPLES <- 8       # eligible samples needed for a coupling test / GSVA
MIN_PAIRS   <- 4       # patients with both lesional and non-lesional eligible
MAX_UNASSIGNED <- 0.40

zs <- function(v) as.numeric(scale(v))

log_cpm <- function(counts) {         # genes x units -> log2(CPM + 1), low genes dropped
  cpm  <- t(t(counts) / colSums(counts)) * 1e6
  keep <- rowSums(cpm > 1) >= 3
  log2(cpm[keep, , drop = FALSE] + 1)
}

panel_score <- function(x, up, down, min_genes = 3) {
  up <- intersect(up, rownames(x)); down <- intersect(down, rownames(x))
  if (length(up) + length(down) < min_genes) return(rep(NA_real_, ncol(x)))
  z <- function(g) t(apply(x[g, , drop = FALSE], 1, zs))
  parts <- NULL
  if (length(up))   parts <- rbind(parts,  z(up))
  if (length(down)) parts <- rbind(parts, -z(down))
  colMeans(parts)
}

classify <- function(beta, p, beta_ref) {
  ifelse(beta > 0 & p < 0.05, "coupled",
         ifelse(p >= 0.05 & beta < 0.5 * beta_ref, "not_coupled", "inconclusive"))
}

# marker-based assignment (plan Section 3): highest mean log1p(CP10K) of a class's markers,
# needs >= 0.5 and a gap >= 0.25 over the runner-up, otherwise "unassigned"
assign_cells <- function(m, markers = MARKERS, thr_min = 0.5, thr_gap = 0.25) {
  lib <- Matrix::colSums(m)
  sc <- vapply(names(markers), function(cl) {
    g <- intersect(markers[[cl]], rownames(m))
    if (length(g) < 3) return(rep(-Inf, ncol(m)))
    x <- as.matrix(m[g, , drop = FALSE])
    colMeans(log1p(sweep(x, 2, lib, "/") * 1e4))
  }, numeric(ncol(m)))
  if (is.null(dim(sc))) sc <- matrix(sc, ncol = length(markers))
  colnames(sc) <- names(markers)
  idx    <- cbind(seq_len(nrow(sc)), max.col(sc, ties.method = "first"))
  top    <- sc[idx]
  sc2    <- sc; sc2[idx] <- -Inf
  second <- apply(sc2, 1, max)
  ifelse(top >= thr_min & (top - second) >= thr_gap, colnames(sc)[idx[, 2]], "unassigned")
}

# read one 10x sample (gz files), collapse rows sharing a gene symbol
read_sample <- function(gsm, name, data_dir) {
  pre  <- file.path(data_dir, paste0(gsm, "_", name))
  m    <- as(readMM(gzfile(paste0(pre, "_matrix.mtx.gz"))), "CsparseMatrix")
  feat <- read.delim(gzfile(paste0(pre, "_features.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)
  bc   <- read.delim(gzfile(paste0(pre, "_barcodes.tsv.gz")), header = FALSE, stringsAsFactors = FALSE)[[1]]
  sym  <- feat[[2]]
  stopifnot(nrow(m) == length(sym), ncol(m) == length(bc))
  f <- factor(sym, levels = sort(unique(sym)))
  m <- fac2sparse(f) %*% m
  rownames(m) <- levels(f); colnames(m) <- bc
  m
}

# plan Section 3: cell filter, then class assignment
process_sample <- function(m) {
  ngene <- Matrix::colSums(m > 0); lib <- Matrix::colSums(m)
  mt <- grep("^MT-", rownames(m))
  mtfrac <- if (length(mt)) Matrix::colSums(m[mt, , drop = FALSE]) / pmax(lib, 1) else rep(0, ncol(m))
  keep <- lib > 0 & ngene >= 200 & mtfrac < 0.20
  m <- m[, keep, drop = FALSE]
  list(m = m, cls = assign_cells(m), n_raw = length(keep), n_kept = sum(keep))
}

# scores for a log2-CPM matrix (units in columns); S2, S3 were dropped in analysis 04
score_units <- function(x) {
  sets <- gs[c("S0", "E2F", "IFNC")]
  n_present <- sapply(sets, function(g) sum(g %in% rownames(x)))
  ok  <- n_present >= 10
  out <- data.frame(row.names = colnames(x))
  if (ncol(x) >= MIN_SAMPLES && any(ok)) {
    gv <- gsva(gsvaParam(x, sets[ok], kcdf = "Gaussian", minSize = 10))
    for (nm in rownames(gv)) out[[nm]] <- as.numeric(gv[nm, colnames(x)])
  }
  for (nm in c("S0", "E2F", "IFNC")) if (!(nm %in% names(out))) out[[nm]] <- NA_real_
  out$S1 <- panel_score(x, panel_up, panel_dn)
  out
}

# standardised slope of S0 and S1 on IFN_c, with E2F sensitivity and classification
couple <- function(sc, meta, extra, label) {
  d <- cbind(sc, meta)
  if (all(is.na(d$IFNC))) return(NULL)
  d$zIFN <- zs(d$IFNC); d$zE2F <- zs(d$E2F)
  rows <- lapply(c("S0", "S1"), function(s) {
    if (all(is.na(d[[s]]))) return(NULL)
    d$zS <- zs(d[[s]])
    m <- tryCatch(lm(as.formula(paste("zS ~ zIFN", extra)), data = d), error = function(e) NULL)
    if (is.null(m)) return(NULL)
    cf <- summary(m)$coefficients
    if (!("zIFN" %in% rownames(cf))) return(NULL)
    ci <- confint(m)["zIFN", ]
    me <- tryCatch(lm(as.formula(paste("zS ~ zIFN", extra, "+ zE2F")), data = d), error = function(e) NULL)
    cfe <- if (!is.null(me) && "zIFN" %in% rownames(summary(me)$coefficients))
      summary(me)$coefficients["zIFN", ] else rep(NA_real_, 4)
    data.frame(analysis = label, score = s, n = nrow(d), beta = cf["zIFN", 1],
               ci_low = ci[[1]], ci_high = ci[[2]], p = cf["zIFN", 4],
               beta_E2Fadj = cfe[[1]], p_E2Fadj = cfe[[4]])
  })
  r <- do.call(rbind, rows)
  if (is.null(r)) return(NULL)
  ref <- r$beta[r$score == "S0"]
  r$class <- if (length(ref) == 1) classify(r$beta, r$p, ref) else NA_character_
  r
}

grp_term <- function(g) if (nlevels(droplevels(g)) > 1) "+ group" else ""

flag_of <- function(tab, classes, score, what = "coupled") {
  r <- tab[tab$cell_type %in% classes & tab$score == score, ]
  if (!nrow(r)) return("not testable")
  if (any(r$class == what, na.rm = TRUE)) "yes" else "no"
}

main_analysis <- function(st, data_dir, outdir) {
  # ---- Step 0: read, filter, assign, pseudobulk ----
  sums <- list(); ncell <- list(); whole <- list(); comp <- list()
  for (i in seq_len(nrow(st))) {
    cat("Reading", st$name[i], "...\n")
    ps <- process_sample(read_sample(st$gsm[i], st$name[i], data_dir))
    m  <- ps$m
    mat <- matrix(0, nrow(m), length(CLASSES), dimnames = list(rownames(m), CLASSES))
    for (cl in CLASSES) { idx <- which(ps$cls == cl); if (length(idx)) mat[, cl] <- Matrix::rowSums(m[, idx, drop = FALSE]) }
    sums[[st$name[i]]]  <- mat
    whole[[st$name[i]]] <- Matrix::rowSums(m)
    ncell[[st$name[i]]] <- table(factor(ps$cls, levels = CLASSES))
    comp[[i]] <- data.frame(sample = st$name[i], group = st$group[i], cells_raw = ps$n_raw,
                            cells_kept = ps$n_kept, t(as.matrix(ncell[[i]])), check.names = FALSE)
    rm(m, ps)
  }
  comp <- do.call(rbind, comp); rownames(comp) <- NULL
  write.csv(comp, file.path(outdir, "cells_per_sample_and_class.csv"), row.names = FALSE)

  unassigned <- sum(comp$unassigned) / sum(comp$cells_kept)
  cat("\nComposition check (cells per class; no scores computed yet):\n"); print(comp, row.names = FALSE)
  cat("\nTotal cells kept:", sum(comp$cells_kept), "| fraction unassigned:", round(unassigned, 3), "\n")
  tot <- colSums(comp[, CLASSES]); cat("Cells per class, all samples:\n"); print(tot)
  pc <- ggplot(data.frame(sample = rep(comp$sample, length(CLASSES)), class = rep(CLASSES, each = nrow(comp)),
                          frac = unlist(lapply(CLASSES, function(cl) comp[[cl]] / comp$cells_kept))),
               aes(sample, frac, fill = class)) + geom_col() + coord_flip() +
    labs(title = "Cell-class composition by sample", x = NULL, y = "fraction of cells")
  ggsave(file.path(outdir, "composition.png"), pc, width = 7, height = 6, dpi = 200)
  if (unassigned > MAX_UNASSIGNED) {
    cat("\n*** STOP: more than ", MAX_UNASSIGNED * 100, "% of cells are unassigned. Per the plan, thresholds may be",
        " revised once and logged BEFORE any scoring. Nothing was scored. ***\n", sep = "")
    return(invisible(NULL))
  }

  common <- Reduce(intersect, lapply(sums, rownames))
  cat("Genes common to all samples:", length(common), "\n")
  meta_all <- data.frame(group = factor(st$group, levels = c("healthy", "nonlesional", "lesional")),
                         status = factor(ifelse(st$group %in% c("lesional", "nonlesional"), st$group, NA),
                                         levels = c("nonlesional", "lesional")),
                         patient = factor(st$patient), row.names = st$name)
  paired_pat <- intersect(st$patient[st$group == "lesional"], st$patient[st$group == "nonlesional"])
  cat("Patients with both lesional and non-lesional samples:", paste(paired_pat, collapse = ", "), "\n")

  # ---- Part 1: whole-sample replication ----
  cat("\n================ PART 1: whole-sample pseudobulk ================\n")
  wc <- sapply(st$name, function(n) whole[[n]][common]); colnames(wc) <- st$name
  x1 <- log_cpm(wc); sc1 <- score_units(x1)
  write.csv(cbind(meta_all, sc1), file.path(outdir, "part1_scores.csv"))
  r1 <- couple(sc1, meta_all, "+ group", "P1 primary (all 17, adj group)")
  sub <- st$name[st$group %in% c("lesional", "nonlesional") & st$patient %in% paired_pat]
  r1b <- couple(sc1[sub, ], droplevels(meta_all[sub, ]), "+ status + patient", "P1 secondary (paired patients)")
  p1 <- rbind(r1, r1b)
  cat("\nPart 1 results (standardised slope of score on IFN_c):\n"); print(p1, digits = 3, row.names = FALSE)
  write.csv(p1, file.path(outdir, "part1_results.csv"), row.names = FALSE)
  c0 <- r1$class[r1$score == "S0"]; c1 <- r1$class[r1$score == "S1"]
  outcome <- if (length(c1) == 1 && identical(c1, "coupled")) "P1-contradicts"
             else if (length(c0) == 1 && length(c1) == 1 && identical(c0, "coupled") && identical(c1, "not_coupled")) "P1-consistent"
             else "P1-inconclusive"
  cat("\nPart 1 outcome:", outcome, "\n")
  d1 <- cbind(sc1, meta_all)
  pp <- ggplot(rbind(data.frame(d1, score_name = "S0 (SenMayo)", value = d1$S0),
                     data.frame(d1, score_name = "S1 (arrest panel)", value = d1$S1)),
               aes(IFNC, value, color = group)) + geom_point(size = 2.5) +
    facet_wrap(~ score_name, scales = "free_y") +
    labs(title = "Part 1: whole-sample scores vs IFN_c (GSE298871)", x = "IFN_c score", y = "score")
  ggsave(file.path(outdir, "part1_scatter.png"), pp, width = 8, height = 4, dpi = 200)

  # ---- Part 2: cell-type resolution ----
  cat("\n================ PART 2: cell-type pseudobulk ================\n")
  elig_rows <- list(); coup <- list(); pdiff <- list(); hnl <- list()
  for (cl in PRIMARY) {
    n_cl <- sapply(st$name, function(n) as.numeric(ncell[[n]][cl]))
    el   <- st$name[n_cl >= MIN_CELLS]
    g_el <- table(factor(st$group[st$name %in% el], levels = c("healthy", "nonlesional", "lesional")))
    pairs_el <- sum(sapply(paired_pat, function(p)
      all(c("lesional", "nonlesional") %in% st$group[st$patient == p & st$name %in% el])))
    elig_rows[[cl]] <- data.frame(cell_type = cl, eligible_samples = length(el),
                                  healthy = g_el[["healthy"]], nonlesional = g_el[["nonlesional"]],
                                  lesional = g_el[["lesional"]], eligible_pairs = pairs_el,
                                  coupling_testable = length(el) >= MIN_SAMPLES,
                                  paired_testable = pairs_el >= MIN_PAIRS)
    cat("\n--", cl, ": eligible samples =", length(el), "| pairs =", pairs_el, "\n")
    if (length(el) < MIN_SAMPLES) { cat("   not testable (fewer than", MIN_SAMPLES, "eligible samples)\n"); next }
    cnt <- sapply(el, function(n) sums[[n]][common, cl]); colnames(cnt) <- el
    x <- log_cpm(cnt); sc <- score_units(x)
    meta <- droplevels(meta_all[el, ])
    r <- couple(sc, meta, grp_term(meta$group), paste("P2", cl))
    if (!is.null(r)) { r$cell_type <- cl; coup[[cl]] <- r }
    if (pairs_el >= MIN_PAIRS) {
      pats <- paired_pat[sapply(paired_pat, function(p)
        all(c("lesional", "nonlesional") %in% st$group[st$patient == p & st$name %in% el]))]
      ss <- el[st$patient[match(el, st$name)] %in% pats & st$group[match(el, st$name)] %in% c("lesional", "nonlesional")]
      dd <- cbind(sc[ss, ], droplevels(meta_all[ss, ]))
      pdiff[[cl]] <- do.call(rbind, lapply(c("IFNC", "S0", "S1"), function(s) {
        if (all(is.na(dd[[s]]))) return(NULL)
        cf <- summary(lm(dd[[s]] ~ status + patient, data = dd))$coefficients
        data.frame(cell_type = cl, score = s, n_pairs = pairs_el,
                   lesional_minus_nonlesional = cf["statuslesional", 1], p = cf["statuslesional", 4])
      }))
    }
    if (g_el[["healthy"]] >= 3 && g_el[["nonlesional"]] >= 3) {
      hs <- el[st$group[match(el, st$name)] %in% c("healthy", "nonlesional")]
      dh <- cbind(sc[hs, ], droplevels(meta_all[hs, ])); dh$group <- droplevels(dh$group)
      hnl[[cl]] <- do.call(rbind, lapply(c("IFNC", "S0", "S1"), function(s) {
        if (all(is.na(dh[[s]]))) return(NULL)
        cf <- summary(lm(dh[[s]] ~ group, data = dh))$coefficients
        data.frame(cell_type = cl, score = s, nonlesional_minus_healthy = cf["groupnonlesional", 1], p = cf["groupnonlesional", 4])
      }))
    }
  }
  elig <- do.call(rbind, elig_rows); rownames(elig) <- NULL
  cat("\nEligibility (primary cell types):\n"); print(elig, row.names = FALSE)
  write.csv(elig, file.path(outdir, "part2_eligibility.csv"), row.names = FALSE)

  tab <- if (length(coup)) do.call(rbind, coup) else NULL
  if (!is.null(tab)) {
    rownames(tab) <- NULL; tab$p_BH <- p.adjust(tab$p, method = "BH")
    cat("\nPart 2 coupling (standardised slope of score on IFN_c within cell type):\n")
    print(tab[, c("cell_type", "score", "n", "beta", "ci_low", "ci_high", "p", "p_BH", "class", "beta_E2Fadj", "p_E2Fadj")],
          digits = 3, row.names = FALSE)
    write.csv(tab, file.path(outdir, "part2_coupling.csv"), row.names = FALSE)
    pf <- ggplot(tab, aes(beta, paste(cell_type, score), color = score)) + geom_vline(xintercept = 0, linetype = 2) +
      geom_pointrange(aes(xmin = ci_low, xmax = ci_high)) +
      labs(title = "Part 2: standardised slope on IFN_c by cell type", x = "slope (95% CI)", y = NULL)
    ggsave(file.path(outdir, "part2_forest.png"), pf, width = 7, height = 4, dpi = 200)
  } else tab <- data.frame(cell_type = character(), score = character(), class = character())

  f1 <- flag_of(tab, c("T_cell", "myeloid"), "S0")
  f2 <- flag_of(tab, c("keratinocyte", "fibroblast"), "S0")
  f3 <- flag_of(tab, PRIMARY, "S1")
  cat("\nFlags:\n F1 (S0 coupled in T cells or myeloid):", f1, "\n F2 (S0 coupled in keratinocytes or fibroblasts):", f2,
      "\n F3 (S1 coupled in any eligible primary class):", f3, "\n")
  interp <- if (f1 == "yes" && f2 == "no") "F1 yes, F2 no: consistent with the coupling reflecting inflammatory-cell content."
            else if (f2 == "yes") "F2 yes: the coupling is not confined to immune cells."
            else "Neither F1-without-F2 nor F2: inconclusive."
  cat("Interpretation (fixed in the plan):", interp, "\n")
  if (f3 == "yes") cat("F3 is yes: notable; arrest genes can follow proliferation or stress changes (see E2F-adjusted columns).\n")
  write.csv(data.frame(flag = c("F1", "F2", "F3", "Part1"), value = c(f1, f2, f3, outcome)),
            file.path(outdir, "flags_and_outcomes.csv"), row.names = FALSE)

  if (length(pdiff)) { pdf_ <- do.call(rbind, pdiff); rownames(pdf_) <- NULL
    cat("\nSecondary: lesional minus non-lesional, paired by patient (no verdict):\n"); print(pdf_, digits = 3, row.names = FALSE)
    write.csv(pdf_, file.path(outdir, "part2_paired_diff.csv"), row.names = FALSE) }
  if (length(hnl)) { h_ <- do.call(rbind, hnl); rownames(h_) <- NULL
    cat("\nSecondary: non-lesional minus healthy (no verdict):\n"); print(h_, digits = 3, row.names = FALSE)
    write.csv(h_, file.path(outdir, "part2_healthy_vs_nonlesional.csv"), row.names = FALSE) }
  invisible(list(comp = comp, part1 = p1, part2 = tab))
}

## ==== END FUNCTIONS ====

## ---------------- sample table from the GEO annotation ----------------
# Sample titles as listed on the GEO page of GSE298871 (copied from the GEO listing). They are used if the
# live GEO lookup fails, which happens occasionally; when the lookup works, it is compared against them.
emb <- data.frame(
  gsm = paste0("GSM", 9025365:9025381),
  title = c("465 (healthy skin)", "486_LS (lesional vitiligo)", "486_NL (non-lesional vitiligo)",
            "632 (healthy skin)", "658 (healthy skin)", "BIO023_LS (lesional vitiligo)",
            "CON161 (healthy skin)", "CON171 (healthy skin)", "CON194 (healthy skin)",
            "DER114_LS (lesional vitiligo)", "DER114_NL (non-lesional vitiligo)",
            "DER115_LS (lesional vitiligo)", "DER115_NL (non-lesional vitiligo)",
            "DER119_LS (lesional vitiligo)", "DER119_NL (non-lesional vitiligo)",
            "DER120_LS (lesional vitiligo)", "DER120_NL (non-lesional vitiligo)"),
  stringsAsFactors = FALSE)
live <- tryCatch({
  ph <- pData(getGEO("GSE298871", GSEMatrix = TRUE)[[1]])
  data.frame(gsm = rownames(ph), title = as.character(ph$title), stringsAsFactors = FALSE)
}, error = function(err) NULL)
if (is.null(live)) {
  cat("GEO lookup failed; using the sample titles embedded in this script.\n")
  tt <- emb
} else {
  tt <- live
  if (!identical(sort(paste(tt$gsm, tt$title)), sort(paste(emb$gsm, emb$title))))
    warning("Live GEO titles differ from the embedded list; check the sample table below.")
}
ttl <- tt$title
st <- data.frame(
  gsm   = tt$gsm,
  name  = sub("^([^ ]+).*$", "\\1", ttl),
  group = ifelse(grepl("non-lesional", ttl, ignore.case = TRUE), "nonlesional",
          ifelse(grepl("lesional", ttl, ignore.case = TRUE), "lesional",
          ifelse(grepl("healthy", ttl, ignore.case = TRUE), "healthy", NA))),
  stringsAsFactors = FALSE)
st$patient <- sub("_(LS|NL)$", "", st$name)
cat("Samples by group (expected 6 healthy, 6 lesional, 5 non-lesional):\n"); print(table(st$group, useNA = "ifany"))
stopifnot(!anyNA(st$group), all(table(st$group)[c("healthy", "lesional", "nonlesional")] == c(6, 6, 5)))
files <- file.path(DATA_DIR, paste0(st$gsm, "_", st$name, c("_matrix.mtx.gz")))
stopifnot(all(file.exists(files)))

## ---------------- frozen scores (same definitions as script 04) ----------------
h  <- msigdbr(species = "Homo sapiens", collection = "H")
c2 <- msigdbr(species = "Homo sapiens", collection = "C2")
allsets <- rbind(h[, c("gs_name", "gene_symbol")], c2[, c("gs_name", "gene_symbol")])
S <- split(allsets$gene_symbol, allsets$gs_name)
ifn <- S$HALLMARK_INTERFERON_GAMMA_RESPONSE
sen_old <- list(S$SAUL_SEN_MAYO, S$FRIDMAN_SENESCENCE_UP,
                S$REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP)
gs <- list(S0   = setdiff(S$SAUL_SEN_MAYO, ifn),
           E2F  = S$HALLMARK_E2F_TARGETS,
           IFNC = setdiff(ifn, unique(unlist(sen_old))))
panel_up <- setdiff(c("CDKN1A", "CDKN2A", "CDKN2B", "GADD45A"), ifn)
panel_dn <- setdiff("LMNB1", ifn)
cat("Set sizes (analysis 04: S0 114, E2F 200, IFNC 181):\n"); print(sapply(gs, length))
cat("Panel genes: up =", paste(panel_up, collapse = ", "), "| down =", paste(panel_dn, collapse = ", "), "\n")
stopifnot(length(gs$S0) == 114, length(gs$E2F) == 200, length(gs$IFNC) == 181)

res <- main_analysis(st, DATA_DIR, OUTDIR)

writeLines(capture.output(sessionInfo()), file.path(OUTDIR, "sessionInfo.txt"))
cat("\nDone. Results saved in", OUTDIR, "\n")
sink()
