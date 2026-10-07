# =====================================================================
# Extension of the systematic replication analysis:
#   (1) second replication cohort GSE17536 (not used for gene selection in this design)
#       + fixed-effect meta-analysis of the two replication cohorts (TCGA, GSE17536)
#   (2) MSI-adjusted sensitivity analysis in TCGA
#   (3) sensitivity of the stromal-vs-epithelial comparison to the class definition
#
# NOTE: the PRIMARY pre-specified analysis was the TCGA-only test already run
#       (systematic_replication.R). Everything here is a robustness / extension
#       analysis and must be reported as such.
#
# Needs in the folder: candidates_discovery.csv, compartment_all_genes.csv,
#   GSE17536_series_matrix.txt.gz, TCGA.COADREAD.sampleMap_HiSeqV2.gz,
#   TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix
# Outputs: systematic2_results.csv, systematic2_summary.txt, systematic2.pdf
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("GEOquery", "Biobase", "AnnotationDbi", "hgu133plus2.db"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
suppressMessages({
  library(GEOquery); library(Biobase); library(AnnotationDbi)
  library(hgu133plus2.db); library(survival)
})

## ---- helpers (copied from dcbld2_analysis.R) --------------------------------------
load_gse <- function(acc) {
  f <- paste0(acc, "_series_matrix.txt.gz")
  if (file.exists(f)) {
    message("Using local file: ", f)
    getGEO(filename = f, GSEMatrix = TRUE, getGPL = FALSE)
  } else {
    message("Downloading ", acc, " from GEO ...")
    getGEO(acc, GSEMatrix = TRUE, getGPL = FALSE)[[1]]
  }
}

to_gene <- function(gse) {
  ex <- exprs(gse)
  if (max(ex, na.rm = TRUE) > 50) ex <- log2(ex + 1)
  sym <- suppressMessages(mapIds(hgu133plus2.db, keys = rownames(ex),
                                 column = "SYMBOL", keytype = "PROBEID",
                                 multiVals = "first"))
  ok <- !is.na(sym)
  ex <- ex[ok, ]; sym <- unname(sym[ok])
  o <- order(rowMeans(ex), decreasing = TRUE)
  ex <- ex[o, ]; sym <- sym[o]
  keep <- !duplicated(sym)
  ex <- ex[keep, ]; rownames(ex) <- sym[keep]
  ex
}

get_os <- function(pd, label) {
  cn <- colnames(pd)
  tcol <- grep("overall.*(time|month|delay)|(^|[^a-z])os.*(delay|time|month)",
               cn, ignore.case = TRUE, value = TRUE)
  ecol <- grep("overall.*event|(^|[^a-z])os.*event|any cause",
               cn, ignore.case = TRUE, value = TRUE)
  if (!length(tcol) || !length(ecol)) stop("Could not find OS columns in ", label)
  message(label, ": time column  = ", tcol[1])
  message(label, ": event column = ", ecol[1])
  raw <- tolower(trimws(as.character(pd[[ecol[1]]])))
  ev <- ifelse(grepl("^(no death|alive|no|0|censored|false)", raw), 0L,
        ifelse(grepl("^(death|dead|deceased|yes|1|true)", raw), 1L, NA_integer_))
  tm <- suppressWarnings(as.numeric(as.character(pd[[tcol[1]]])))
  d <- data.frame(sample = rownames(pd), time = tm, event = ev)
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
  message(label, ": patients = ", nrow(d), ", events = ", sum(d$event))
  d
}

clean <- function(x) {
  x <- trimws(as.character(x))
  x[tolower(x) %in% c("n/a", "na", "", "nan", "unknown", "not available")] <- NA
  x
}
pick <- function(pd, pattern) {
  cn <- grep(pattern, colnames(pd), ignore.case = TRUE, value = TRUE)
  if (length(cn)) cn[1] else NA_character_
}

# age, sex, stage (I-II vs III-IV), MMR status, when the cohort has them
make_cov <- function(pd, samples, label) {
  rn <- match(samples, rownames(pd))
  df <- data.frame(sample = samples)
  a <- pick(pd, "^age");  s <- pick(pd, "^(sex|gender)")
  st <- pick(pd, "stage"); m <- pick(pd, "mmr")
  message(label, " covariate columns -> age: ", a, " | sex: ", s,
          " | stage: ", st, " | MMR: ", m)
  if (!is.na(a))  df$age <- suppressWarnings(as.numeric(clean(pd[[a]])))[rn]
  if (!is.na(s))  df$sex <- factor(tolower(clean(pd[[s]]))[rn])
  if (!is.na(st)) {
    x <- clean(pd[[st]])[rn]
    mm <- regexpr("[0-4]", x)
    num <- ifelse(!is.na(mm) & mm > 0, as.integer(substr(x, mm, mm)), NA_integer_)
    df$stage <- factor(ifelse(num <= 2, "I-II", "III-IV"), levels = c("I-II", "III-IV"))
  }
  if (!is.na(m))  df$mmr <- factor(clean(pd[[m]])[rn])
  df
}


find1 <- function(pat) {
  f <- list.files(pattern = pat, ignore.case = TRUE)
  f <- f[!grepl("\\.(R|csv|pdf|png|txt)$", f)]
  if (!length(f)) stop("File matching '", pat, "' not found in this folder.")
  f[1]
}
# estimates: c(logHR_uni, se_uni, p_uni, logHR_adj, se_adj, p_adj) for one gene
cox_est <- function(v, dd, covs) {
  out <- rep(NA_real_, 6)
  if (is.na(sd(v)) || sd(v) == 0) return(out)
  dd$z <- as.numeric(scale(v))
  r1 <- tryCatch(summary(coxph(Surv(time, event) ~ z, data = dd)), error = function(e) NULL)
  if (!is.null(r1)) out[1:3] <- c(r1$coefficients[, "coef"], r1$coefficients[, "se(coef)"],
                                  r1$coefficients[, "Pr(>|z|)"])
  if (length(covs)) {
    f <- as.formula(paste("Surv(time, event) ~ z +", paste(covs, collapse = " + ")))
    r2 <- tryCatch(summary(coxph(f, data = dd)), error = function(e) NULL)
    if (!is.null(r2)) out[4:6] <- c(r2$coefficients["z", "coef"], r2$coefficients["z", "se(coef)"],
                                    r2$coefficients["z", "Pr(>|z|)"])
  }
  out
}
informative <- function(dd, covs)
  covs[vapply(covs, function(k) mean(is.na(dd[[k]])) < 0.3 && length(unique(na.omit(dd[[k]]))) > 1, logical(1))]

disc <- read.csv("candidates_discovery.csv", stringsAsFactors = FALSE)
comp <- read.csv("compartment_all_genes.csv", stringsAsFactors = FALSE)
disc <- disc[!is.na(disc$FDR_disc) & disc$FDR_disc < 0.05, ]
message("Discovery genes with FDR < 0.05: ", nrow(disc))

## ---- TCGA --------------------------------------------------------------------------
f_expr <- find1("HiSeqV2"); f_clin <- find1("clinicalMatrix")
ex <- read.delim(gzfile(f_expr), check.names = FALSE, stringsAsFactors = FALSE)
rownames(ex) <- ex[[1]]; ex <- ex[, -1, drop = FALSE]
tum <- colnames(ex)[substr(colnames(ex), 14, 15) == "01"]
cl <- read.delim(f_clin, check.names = FALSE, stringsAsFactors = FALSE)
g1 <- function(p) { k <- grep(p, colnames(cl), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
cv <- g1("^vital_status$"); cd <- g1("^days_to_death$"); cf <- g1("^days_to_last_followup$")
dead <- toupper(trimws(cl[[cv]])) == "DECEASED"
tdead <- suppressWarnings(as.numeric(cl[[cd]])); tlast <- suppressWarnings(as.numeric(cl[[cf]]))
d <- data.frame(sample = cl[[1]], time = ifelse(dead, tdead, tlast) / 30.4375, event = as.integer(dead),
                stringsAsFactors = FALSE)
a <- g1("^age_at_initial"); gg <- g1("^gender$"); s <- g1("^pathologic_stage$"); mi <- g1("^microsatellite_instability$")
if (!is.na(a))  d$age <- suppressWarnings(as.numeric(cl[[a]]))
if (!is.na(gg)) d$gender <- factor(cl[[gg]])
if (!is.na(s))  d$stage <- factor(ifelse(grepl("Stage (III|IV)", cl[[s]]), "III-IV",
                                  ifelse(grepl("Stage (I|II)", cl[[s]]), "I-II", NA)), levels = c("I-II", "III-IV"))
if (!is.na(mi)) { m0 <- toupper(trimws(cl[[mi]])); m0[m0 %in% c("", "NA", "[NOT AVAILABLE]", "[UNKNOWN]")] <- NA
                  d$msi <- factor(m0) }
d <- d[d$sample %in% tum & !is.na(d$time) & d$time > 0 & !is.na(d$event), ]
covs_t <- informative(d, intersect(c("age", "gender", "stage"), colnames(d)))
msi_n <- if ("msi" %in% colnames(d)) sum(!is.na(d$msi)) else 0
covs_m <- if (msi_n >= 100 && length(unique(na.omit(d$msi))) > 1) c(covs_t, "msi") else NULL
message("TCGA patients: ", nrow(d), ", deaths: ", sum(d$event), " | covariates: ", paste(covs_t, collapse = ", "),
        " | patients with MSI status: ", msi_n)
X <- as.matrix(ex[, d$sample, drop = FALSE])
genes_t <- intersect(disc$gene, rownames(X))
mt <- t(vapply(genes_t, function(g) cox_est(X[g, ], d, covs_t), numeric(6)))
colnames(mt) <- c("b_t", "se_t", "p_t", "b_ta", "se_ta", "p_ta")
rt <- data.frame(gene = genes_t, mt, row.names = NULL)
if (!is.null(covs_m)) {
  mm <- t(vapply(genes_t, function(g) cox_est(X[g, ], d, covs_m)[4:6], numeric(3)))
  colnames(mm) <- c("b_tm", "se_tm", "p_tm"); rt <- cbind(rt, mm)
} else { rt$b_tm <- NA_real_; rt$se_tm <- NA_real_; rt$p_tm <- NA_real_ }

## ---- GSE17536 ---------------------------------------------------------------------------
g2 <- load_gse("GSE17536"); pd2 <- pData(g2); ex2 <- to_gene(g2)
os2 <- get_os(pd2, "GSE17536"); os2 <- os2[os2$sample %in% colnames(ex2), ]
cv2 <- make_cov(pd2, os2$sample, "GSE17536")
d2 <- if (ncol(cv2) > 1) cbind(os2, cv2[, -1, drop = FALSE]) else os2
covs_g <- informative(d2, setdiff(colnames(cv2), "sample"))
message("GSE17536 patients: ", nrow(d2), ", deaths: ", sum(d2$event), " | covariates: ", paste(covs_g, collapse = ", "))
genes_g <- intersect(disc$gene, rownames(ex2))
mg <- t(vapply(genes_g, function(g) cox_est(ex2[g, d2$sample], d2, covs_g), numeric(6)))
colnames(mg) <- c("b_g", "se_g", "p_g", "b_ga", "se_ga", "p_ga")
rg <- data.frame(gene = genes_g, mg, row.names = NULL)

## ---- merge, meta-analysis, flags ----------------------------------------------------------
res <- merge(disc[, c("gene", "HR_disc", "FDR_disc")], rt, by = "gene")
res <- merge(res, rg, by = "gene")
res <- merge(res, comp, by = "gene")
res <- res[!is.na(res$b_t) & !is.na(res$b_g), ]
meta2 <- function(b1, s1, b2, s2) {
  w1 <- 1 / s1^2; w2 <- 1 / s2^2
  b <- (w1 * b1 + w2 * b2) / (w1 + w2); se <- sqrt(1 / (w1 + w2))
  list(b = b, p = 2 * pnorm(-abs(b / se)))
}
mu <- meta2(res$b_t, res$se_t, res$b_g, res$se_g); res$b_mu <- mu$b; res$p_mu <- mu$p
ma <- meta2(res$b_ta, res$se_ta, res$b_ga, res$se_ga); res$b_ma <- ma$b; res$p_ma <- ma$p
sg <- ifelse(res$HR_disc > 1, 1, -1)
okdir <- function(b) sign(b) == sg
res$rep_tcga_uni  <- !is.na(res$p_t)  & res$p_t  < 0.05 & okdir(res$b_t)
res$rep_tcga_adj  <- !is.na(res$p_ta) & res$p_ta < 0.05 & okdir(res$b_ta)
res$rep_gse_uni   <- !is.na(res$p_g)  & res$p_g  < 0.05 & okdir(res$b_g)
res$rep_gse_adj   <- !is.na(res$p_ga) & res$p_ga < 0.05 & okdir(res$b_ga)
res$rep_meta_uni  <- !is.na(res$p_mu) & res$p_mu < 0.05 & okdir(res$b_mu)
res$rep_meta_adj  <- !is.na(res$p_ma) & res$p_ma < 0.05 & okdir(res$b_ma)
res$rep_tcga_msi  <- !is.na(res$p_tm) & res$p_tm < 0.05 & okdir(res$b_tm)
write.csv(res, "systematic2_results.csv", row.names = FALSE)

## ---- summary --------------------------------------------------------------------------------
classify <- function(df, thr, prefix) {
  ifelse(df[[paste0(prefix, "stromal")]] > thr, "stromal",
  ifelse(df[[paste0(prefix, "epithelial")]] > thr, "epithelial",
  ifelse(df[[paste0(prefix, "immune")]] > thr, "immune", "mixed")))
}
test_class <- function(df, thr, prefix, outcome) {
  k <- classify(df, thr, prefix)
  s <- df[[outcome]][k == "stromal"]; e <- df[[outcome]][k == "epithelial"]
  p <- tryCatch(fisher.test(matrix(c(sum(e), length(e) - sum(e), sum(s), length(s) - sum(s)), 2))$p.value,
                error = function(e) NA)
  data.frame(basis = prefix, threshold = thr, outcome = outcome, n_stromal = length(s), rep_stromal = sum(s),
             n_epithelial = length(e), rep_epithelial = sum(e), fisher_p = round(p, 4))
}
sink("systematic2_summary.txt", split = TRUE)
cat("Genes analysed (present in TCGA, GSE17536 and the single-cell table):", nrow(res), "\n\n")
rates <- sapply(c("rep_tcga_uni", "rep_tcga_adj", "rep_gse_uni", "rep_gse_adj", "rep_meta_uni", "rep_meta_adj",
                  "rep_tcga_msi"), function(k) c(n = sum(res[[k]]), pct = round(100 * mean(res[[k]]), 1)))
print(t(rates))
cat("\nMedian risk-aligned log(HR): discovery", round(median(sg * log(res$HR_disc)), 3),
    "| TCGA", round(median(sg * res$b_t), 3), "| GSE17536", round(median(sg * res$b_g), 3),
    "| pooled", round(median(sg * res$b_mu), 3), "\n")
cat("Same direction as discovery: TCGA", sum(okdir(res$b_t)), "| GSE17536", sum(okdir(res$b_g)),
    "| both", sum(okdir(res$b_t) & okdir(res$b_g)), "of", nrow(res), "\n\n")

cat("--- stromal- vs epithelial-dominant genes: sensitivity to the class definition ---\n")
tests <- do.call(rbind, lapply(c("rep_tcga_adj", "rep_meta_uni", "rep_meta_adj"), function(o)
  do.call(rbind, c(lapply(c(0.4, 0.5, 0.6), function(t) test_class(res, t, "share_bulk_", o)),
                   list(test_class(res, 0.5, "share_spec_", o))))))
print(tests, row.names = FALSE)

cat("\n--- genes replicated in the pooled adjusted analysis, by class (>50% bulk share) ---\n")
res$class50 <- classify(res, 0.5, "share_bulk_")
print(table(res$class50, res$rep_meta_adj))
cat("\n--- MSI sensitivity: genes replicated after adjustment in TCGA ---\n")
sub <- res[res$rep_tcga_adj, c("gene", "class50", "p_ta", "p_tm", "rep_tcga_msi")]
sub$HR_adj <- round(exp(res$b_ta[res$rep_tcga_adj]), 3); sub$HR_adjMSI <- round(exp(res$b_tm[res$rep_tcga_adj]), 3)
print(sub, row.names = FALSE)
cat("\nGenes replicated in the pooled adjusted analysis:\n")
print(res[res$rep_meta_adj, c("gene", "class50", "HR_disc")], row.names = FALSE)
sink()

pdf("systematic2.pdf", width = 7.5, height = 5.5)
cls <- c("stromal", "epithelial")
pc <- sapply(cls, function(k) sapply(c("rep_meta_uni", "rep_meta_adj"),
                                     function(o) 100 * mean(res[[o]][res$class50 == k])))
nn <- sapply(cls, function(k) sum(res$class50 == k))
barplot(pc, beside = TRUE, names.arg = paste0(cls, "-dominant\n(n=", nn, ")"),
        col = c("grey70", "firebrick"), ylab = "% genes replicated (TCGA + GSE17536 pooled)",
        main = "Replication by cellular compartment (extension analysis)")
legend("topright", c("univariable", "adjusted"), fill = c("grey70", "firebrick"), bty = "n")
dev.off()
cat("Done. See systematic2_summary.txt, systematic2_results.csv, systematic2.pdf\n")
