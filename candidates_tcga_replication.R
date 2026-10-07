# =====================================================================
# Replication of the 18 GEO-derived candidate genes in TCGA-COADREAD (OS)
# Needs in the working folder (all already there):
#   candidates_ranked.csv                       (from crc_geo_candidates.R)
#   TCGA.COADREAD.sampleMap_HiSeqV2.gz          (UCSC Xena expression)
#   TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix
# Outputs: candidates_tcga_replication.csv, candidates_tcga_replication.pdf
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
suppressMessages(library(survival))

find1 <- function(pat) {
  f <- list.files(pattern = pat, ignore.case = TRUE)
  f <- f[!grepl("\\.(R|csv|pdf|png)$", f)]
  if (!length(f)) stop("File matching '", pat, "' not found in this folder.")
  f[1]
}
cand   <- read.csv("candidates_ranked.csv", stringsAsFactors = FALSE)
f_expr <- find1("HiSeqV2"); f_clin <- find1("clinicalMatrix")

## ---- expression (tumors only) -------------------------------------------------
ex <- read.delim(gzfile(f_expr), check.names = FALSE, stringsAsFactors = FALSE)
rownames(ex) <- ex[[1]]; ex <- ex[, -1, drop = FALSE]
tum <- colnames(ex)[substr(colnames(ex), 14, 15) == "01"]

## ---- overall survival + covariates from the clinicalMatrix ---------------------------
cl <- read.delim(f_clin, check.names = FALSE, stringsAsFactors = FALSE)
g1 <- function(p) { k <- grep(p, colnames(cl), ignore.case = TRUE, value = TRUE); if (length(k)) k[1] else NA }
cv <- g1("^vital_status$"); cd <- g1("^days_to_death$"); cf <- g1("^days_to_last_followup$")
if (anyNA(c(cv, cd, cf))) stop("vital_status / days_to_death / days_to_last_followup not found")
dead <- toupper(trimws(cl[[cv]])) == "DECEASED"
tdead <- suppressWarnings(as.numeric(cl[[cd]])); tlast <- suppressWarnings(as.numeric(cl[[cf]]))
d <- data.frame(sample = cl[[1]], time = ifelse(dead, tdead, tlast) / 30.4375, event = as.integer(dead),
                stringsAsFactors = FALSE)
a <- g1("^age_at_initial"); gg <- g1("^gender$"); s <- g1("^pathologic_stage$")
if (!is.na(a))  d$age <- suppressWarnings(as.numeric(cl[[a]]))
if (!is.na(gg)) d$gender <- factor(cl[[gg]])
if (!is.na(s))  d$stage <- factor(ifelse(grepl("Stage (III|IV)", cl[[s]]), "III-IV",
                                  ifelse(grepl("Stage (I|II)", cl[[s]]), "I-II", NA)),
                                  levels = c("I-II", "III-IV"))
d <- d[d$sample %in% tum & !is.na(d$time) & d$time > 0 & !is.na(d$event), ]
covs <- intersect(c("age", "gender", "stage"), colnames(d))
covs <- covs[vapply(covs, function(k) mean(is.na(d[[k]])) < 0.3 && length(unique(na.omit(d[[k]]))) > 1,
                    logical(1))]
message("TCGA patients with OS: ", nrow(d), ", deaths: ", sum(d$event), " | covariates: ", paste(covs, collapse = ", "))

## ---- Cox for every candidate ----------------------------------------------------------
one <- function(g) {
  out <- data.frame(gene = g, in_tcga = g %in% rownames(ex), HR_tcga = NA_real_, lo_tcga = NA_real_,
                    hi_tcga = NA_real_, p_tcga = NA_real_, HR_adj = NA_real_, p_adj = NA_real_)
  if (!out$in_tcga) return(out)
  dd <- d
  dd$z <- as.numeric(scale(as.numeric(unlist(ex[g, dd$sample]))))
  r <- tryCatch(summary(coxph(Surv(time, event) ~ z, data = dd)), error = function(e) NULL)
  if (is.null(r)) return(out)
  out$HR_tcga <- r$coefficients[, "exp(coef)"]; out$p_tcga <- r$coefficients[, "Pr(>|z|)"]
  out$lo_tcga <- r$conf.int[, "lower .95"]; out$hi_tcga <- r$conf.int[, "upper .95"]
  if (length(covs)) {
    f <- as.formula(paste("Surv(time, event) ~ z +", paste(covs, collapse = " + ")))
    r2 <- tryCatch(summary(coxph(f, data = dd)), error = function(e) NULL)
    if (!is.null(r2)) { out$HR_adj <- r2$coefficients["z", "exp(coef)"]; out$p_adj <- r2$coefficients["z", "Pr(>|z|)"] }
  }
  out
}
res <- do.call(rbind, lapply(cand$gene, one))
res <- merge(cand[, c("gene", "HR_disc", "HR_val", "pubmed_CRC_hits")], res, by = "gene", sort = FALSE)
res$same_direction <- (res$HR_tcga > 1) == (res$HR_disc > 1)
res$FDR_tcga <- p.adjust(res$p_tcga, "BH")
res$replicated_uni <- !is.na(res$p_tcga) & res$p_tcga < 0.05 & res$same_direction
res$replicated_adj <- !is.na(res$p_adj) & res$p_adj < 0.05 & ((res$HR_adj > 1) == (res$HR_disc > 1))
res <- res[order(res$p_tcga), ]
write.csv(res, "candidates_tcga_replication.csv", row.names = FALSE)
print(res[, c("gene", "HR_disc", "HR_val", "HR_tcga", "p_tcga", "HR_adj", "p_adj", "replicated_uni", "replicated_adj")], digits = 3)

tested <- sum(res$in_tcga & !is.na(res$p_tcga))
cat("\nGenes present and tested in TCGA:", tested, "of", nrow(res), "\n")
cat("Replicated (p<0.05, same direction), univariable:", sum(res$replicated_uni), "\n")
cat("Replicated after adjustment (age, gender, stage):", sum(res$replicated_adj), "\n")
sg <- ifelse(res$HR_disc > 1, 1, -1)
ok <- !is.na(res$HR_tcga)
cat(sprintf("Median risk-aligned log(HR): discovery %.3f | second GEO cohort %.3f | TCGA %.3f\n",
            median((sg * log(res$HR_disc))[ok]), median((sg * log(res$HR_val))[ok]),
            median((sg * log(res$HR_tcga))[ok])))

## ---- figure: discovery vs TCGA effect sizes -------------------------------------------------
pdf("candidates_tcga_replication.pdf", width = 7, height = 6.5)
x <- log(res$HR_disc[ok]); y <- log(res$HR_tcga[ok])
plot(x, y, pch = 16, col = ifelse(res$replicated_uni[ok], "firebrick", "grey50"),
     xlab = "log(HR) in discovery cohort (GSE39582)", ylab = "log(HR) in TCGA-COADREAD",
     main = "Replication of GEO-derived prognostic candidates")
abline(h = 0, v = 0, lty = 3); abline(0, 1, lty = 2)
text(x, y, labels = res$gene[ok], pos = 3, cex = 0.7)
legend("topleft", c("replicated (p<0.05)", "not replicated"), pch = 16,
       col = c("firebrick", "grey50"), bty = "n")
dev.off()
cat("Done. See candidates_tcga_replication.csv and candidates_tcga_replication.pdf\n")
