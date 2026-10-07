# =====================================================================
# DCBLD2: independent validation in TCGA-COADREAD (overall survival, OS)
# and an OS-only meta-analysis with GSE39582 and GSE17536.
#
# Put these files (downloaded from UCSC Xena, TCGA Colon and Rectal Cancer,
# COADREAD) in the working folder; the script finds them by name pattern:
#   *HiSeqV2*           gene expression (log2(RSEM+1)), rows = gene symbols
#   *survival*          (optional) survival table (columns: sample, OS, OS.time ...);
#                       if absent, OS is built from clinicalMatrix
#   *clinicalMatrix*    (optional) age, gender, stage for the multivariable model
# Also uses DCBLD2_cox_results.csv (from dcbld2_analysis.R) for the GEO cohorts.
#
# Outputs: DCBLD2_tcga_cox.csv, DCBLD2_meta_OS.csv,
#          DCBLD2_tcga_figures.pdf, DCBLD2_forest_OS.pdf
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
GENE <- "DCBLD2"
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
suppressMessages(library(survival))

find1 <- function(pat, required = TRUE) {
  f <- list.files(pattern = pat, ignore.case = TRUE)
  f <- f[!grepl("\\.(R|csv|pdf|png)$", f)]
  if (!length(f)) {
    if (required) stop("File matching '", pat, "' not found in this folder.")
    return(NA_character_)
  }
  f[1]
}
f_expr <- find1("HiSeqV2"); f_surv <- find1("survival", FALSE); f_clin <- find1("clinicalMatrix", FALSE)
message("Using: ", f_expr, " | ", f_surv, " | ", f_clin)

## ---- expression ------------------------------------------------------------
ex <- read.delim(gzfile(f_expr), check.names = FALSE, stringsAsFactors = FALSE)
rownames(ex) <- ex[[1]]; ex <- ex[, -1, drop = FALSE]
if (!GENE %in% rownames(ex)) stop(GENE, " not found in the expression file")
x <- as.numeric(unlist(ex[GENE, ])); names(x) <- colnames(ex)
stype <- substr(names(x), 14, 15)              # 01 = primary tumor, 11 = normal
message("Samples: tumor(01) = ", sum(stype == "01"), ", normal(11) = ", sum(stype == "11"))

## ---- survival ------------------------------------------------------------------
if (!is.na(f_surv)) {
  sv <- read.delim(gzfile(f_surv), check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample", "OS", "OS.time") %in% colnames(sv)))
    stop("Survival file lacks sample/OS/OS.time. Columns found: ", paste(colnames(sv), collapse = ", "))
} else {
  # no survival file: build overall survival from the clinicalMatrix
  if (is.na(f_clin)) stop("Need either the survival file or the clinicalMatrix file.")
  c0 <- read.delim(f_clin, check.names = FALSE, stringsAsFactors = FALSE)
  g1 <- function(p) { k <- grep(p, colnames(c0), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
  cv <- g1("^vital_status$"); cd <- g1("^days_to_death$"); cf <- g1("^days_to_last_followup$")
  if (anyNA(c(cv, cd, cf))) {
    message("clinicalMatrix columns:\n", paste(colnames(c0), collapse = "\n"))
    stop("Could not find vital_status / days_to_death / days_to_last_followup. Send me the column list above.")
  }
  message("Survival built from clinicalMatrix: ", cv, ", ", cd, ", ", cf)
  dead <- toupper(trimws(c0[[cv]])) == "DECEASED"
  tdead <- suppressWarnings(as.numeric(c0[[cd]])); tlast <- suppressWarnings(as.numeric(c0[[cf]]))
  sv <- data.frame(sample = c0[[1]], OS = as.integer(dead), OS.time = ifelse(dead, tdead, tlast))
}
d <- data.frame(sample = names(x)[stype == "01"], gene = x[stype == "01"], stringsAsFactors = FALSE)
d <- merge(d, sv[, c("sample", "OS", "OS.time")], by = "sample")
d$time  <- d$OS.time / 30.4375                 # days -> months
d$event <- d$OS
d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
message("Patients with OS: ", nrow(d), ", deaths: ", sum(d$event))
d$z <- as.numeric(scale(d$gene))

cox_tab <- function(m, model) {
  s <- summary(m)
  data.frame(model = model, term = rownames(s$coefficients), HR = s$coefficients[, "exp(coef)"],
             lower95 = s$conf.int[, "lower .95"], upper95 = s$conf.int[, "upper .95"],
             p = s$coefficients[, "Pr(>|z|)"], n = m$n, events = m$nevent, row.names = NULL)
}
res <- list(cox_tab(coxph(Surv(time, event) ~ z, data = d), "univariable"))

## optional multivariable model (age, gender, stage)
if (!is.na(f_clin)) {
  cl <- read.delim(f_clin, check.names = FALSE, stringsAsFactors = FALSE)
  idc <- colnames(cl)[1]
  pick <- function(p) { k <- grep(p, colnames(cl), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
  a <- pick("^age_at_initial"); g <- pick("^gender$"); s <- pick("^pathologic_stage$")
  message("Clinical columns: age = ", a, " | gender = ", g, " | stage = ", s)
  cc <- data.frame(sample = cl[[idc]], stringsAsFactors = FALSE)
  if (!is.na(a)) cc$age <- suppressWarnings(as.numeric(cl[[a]]))
  if (!is.na(g)) cc$gender <- factor(cl[[g]])
  if (!is.na(s)) {
    st <- cl[[s]]
    cc$stage <- factor(ifelse(grepl("Stage (III|IV)", st), "III-IV",
                       ifelse(grepl("Stage (I|II)", st), "I-II", NA)), levels = c("I-II", "III-IV"))
  }
  dd <- merge(d, cc, by = "sample")
  covs <- setdiff(colnames(cc), "sample")
  covs <- covs[vapply(covs, function(k) mean(is.na(dd[[k]])) < 0.3 && length(unique(na.omit(dd[[k]]))) > 1, logical(1))]
  if (length(covs)) {
    f <- as.formula(paste("Surv(time, event) ~ z +", paste(covs, collapse = " + ")))
    res[[2]] <- cox_tab(coxph(f, data = dd), paste("multivariable:", paste(covs, collapse = "+")))
  }
}
tcga_res <- do.call(rbind, res)
write.csv(tcga_res, "DCBLD2_tcga_cox.csv", row.names = FALSE)
print(tcga_res[tcga_res$term == "z", ])

## ---- figures -------------------------------------------------------------------------
pdf("DCBLD2_tcga_figures.pdf", width = 7, height = 6)
# tumor vs normal (unpaired and paired)
tn <- data.frame(sample = names(x), v = x, type = stype)
tn <- tn[tn$type %in% c("01", "11"), ]
tn$type <- factor(ifelse(tn$type == "01", "Tumor", "Normal"), levels = c("Normal", "Tumor"))
if (sum(tn$type == "Normal") >= 5) {
  pw <- wilcox.test(v ~ type, data = tn)$p.value
  boxplot(v ~ type, data = tn, col = c("grey80", "salmon"), outline = FALSE,
          ylab = paste(GENE, "log2(RSEM+1)"), main = paste(GENE, "- TCGA-COADREAD"))
  stripchart(v ~ type, data = tn, vertical = TRUE, method = "jitter", pch = 16,
             col = rgb(0, 0, 0, 0.3), add = TRUE)
  legend("topleft", sprintf("Wilcoxon p = %.2g", pw), bty = "n")
}
# Kaplan-Meier, median split
d$grp <- factor(ifelse(d$gene > median(d$gene), "High", "Low"), levels = c("Low", "High"))
fit <- survfit(Surv(time, event) ~ grp, data = d)
lr  <- survdiff(Surv(time, event) ~ grp, data = d)
plot(fit, col = c("steelblue", "firebrick"), lwd = 2, mark.time = TRUE,
     xlab = "Time (months)", ylab = "Overall survival", main = paste(GENE, "- TCGA-COADREAD"))
legend("bottomleft", c("Low", "High"), col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
legend("topright", sprintf("log-rank p = %.2g (median split)", 1 - pchisq(lr$chisq, 1)), bty = "n")
dev.off()

## ---- OS-only meta-analysis ------------------------------------------------------------
meta <- function(b, se) {
  w <- 1 / se^2; k <- length(b)
  bf <- sum(w * b) / sum(w); sef <- sqrt(1 / sum(w))
  Q <- sum(w * (b - bf)^2)
  tau2 <- if (k > 1) max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w))) else 0
  wr <- 1 / (se^2 + tau2); br <- sum(wr * b) / sum(wr); ser <- sqrt(1 / sum(wr))
  I2 <- if (k > 1 && Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0
  est <- c(bf, br); s <- c(sef, ser)
  data.frame(model = c("fixed-effect", "random-effects"), k = k, HR = exp(est),
             lower95 = exp(est - 1.96 * s), upper95 = exp(est + 1.96 * s),
             p = 2 * pnorm(-abs(est / s)), I2_percent = I2)
}
u <- tcga_res[tcga_res$term == "z" & tcga_res$model == "univariable", ]
tc <- data.frame(cohort = "TCGA-COADREAD (independent)", n = u$n, events = u$events,
                 logHR = log(u$HR), se = (log(u$upper95) - log(u$lower95)) / (2 * 1.96), HR = u$HR)
all <- tc
if (file.exists("DCBLD2_cox_results.csv")) {
  cx <- read.csv("DCBLD2_cox_results.csv")
  cx <- cx[cx$term == "gene" & cx$model == "univariable", ]
  geo <- data.frame(cohort = paste0(cx$cohort, " (used for selection)"), n = cx$n, events = cx$events,
                    logHR = log(cx$HR), se = (log(cx$upper95) - log(cx$lower95)) / (2 * 1.96), HR = cx$HR)
  all <- rbind(geo, tc)
}
mt <- cbind(set = "all OS cohorts", meta(all$logHR, all$se))
write.csv(mt, "DCBLD2_meta_OS.csv", row.names = FALSE)
print(mt)

pdf("DCBLD2_forest_OS.pdf", width = 9, height = 4.5)
pooled <- mt[mt$model == "random-effects", ]
fd <- data.frame(label = paste0(all$cohort, " [OS, n=", all$n, "]"), HR = all$HR,
                 lo = exp(all$logHR - 1.96 * all$se), hi = exp(all$logHR + 1.96 * all$se))
fd <- rbind(fd, data.frame(label = "Pooled (random-effects)", HR = pooled$HR,
                           lo = pooled$lower95, hi = pooled$upper95))
k <- nrow(fd); par(mar = c(4, 17, 3, 2))
plot(NA, xlim = c(min(fd$lo, 0.8), max(fd$hi, 1.5)), ylim = c(0.5, k + 0.5), log = "x", yaxt = "n",
     xlab = paste("Hazard ratio per 1 SD of", GENE), ylab = "", main = paste(GENE, "- overall survival"))
abline(v = 1, lty = 2)
segments(fd$lo, k:1, fd$hi, k:1, lwd = 2)
points(fd$HR, k:1, pch = 15)
axis(2, at = k:1, labels = fd$label, las = 1, cex.axis = 0.75)
dev.off()
cat("Done. See DCBLD2_tcga_cox.csv, DCBLD2_meta_OS.csv, DCBLD2_tcga_figures.pdf, DCBLD2_forest_OS.pdf\n")
