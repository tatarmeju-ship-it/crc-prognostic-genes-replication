# =====================================================================
# Does the cellular origin of a bulk prognostic gene predict whether it replicates?
#
# Discovery set : all genes with FDR_disc < 0.05 in GSE39582 (candidates_discovery.csv)
#                 -> selected WITHOUT looking at TCGA or at single-cell data
# Replication   : TCGA-COADREAD overall survival (univariable and adjusted for
#                 age, gender, stage)
# Compartment   : compartment_all_genes.csv (from compartment_all_genes.py)
#
# PRE-SPECIFIED PREDICTION (write it down with a date BEFORE running this script):
#   stromal-dominant genes (>50% of tumor UMIs from stromal cells) replicate after
#   stage adjustment more often than epithelial-dominant genes.
#
# Needs in the folder: candidates_discovery.csv, compartment_all_genes.csv,
#   TCGA.COADREAD.sampleMap_HiSeqV2.gz, TCGA.COADREAD.sampleMap_COADREAD_clinicalMatrix
# Outputs: systematic_replication.csv, systematic_summary.txt, systematic_replication.pdf
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
suppressMessages(library(survival))

find1 <- function(pat) {
  f <- list.files(pattern = pat, ignore.case = TRUE)
  f <- f[!grepl("\\.(R|csv|pdf|png|txt)$", f)]
  if (!length(f)) stop("File matching '", pat, "' not found in this folder.")
  f[1]
}
disc <- read.csv("candidates_discovery.csv", stringsAsFactors = FALSE)
comp <- read.csv("compartment_all_genes.csv", stringsAsFactors = FALSE)
disc <- disc[!is.na(disc$FDR_disc) & disc$FDR_disc < 0.05, ]
message("Discovery genes with FDR < 0.05: ", nrow(disc))

## ---- TCGA expression (tumors) + OS + covariates ---------------------------------------
f_expr <- find1("HiSeqV2"); f_clin <- find1("clinicalMatrix")
ex <- read.delim(gzfile(f_expr), check.names = FALSE, stringsAsFactors = FALSE)
rownames(ex) <- ex[[1]]; ex <- ex[, -1, drop = FALSE]
tum <- colnames(ex)[substr(colnames(ex), 14, 15) == "01"]
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
fadj <- as.formula(paste("Surv(time, event) ~ z +", paste(covs, collapse = " + ")))
message("TCGA patients: ", nrow(d), ", deaths: ", sum(d$event), " | covariates: ", paste(covs, collapse = ", "))
X <- as.matrix(ex[, d$sample, drop = FALSE])

## ---- Cox for every discovery gene -------------------------------------------------------
genes <- intersect(disc$gene, rownames(X))
message("Discovery genes present in TCGA: ", length(genes))
one <- function(g) {
  v <- X[g, ]
  if (is.na(sd(v)) || sd(v) == 0) return(c(NA, NA, NA, NA))
  dd <- d; dd$z <- as.numeric(scale(v))
  r1 <- tryCatch(summary(coxph(Surv(time, event) ~ z, data = dd)), error = function(e) NULL)
  r2 <- tryCatch(summary(coxph(fadj, data = dd)), error = function(e) NULL)
  c(if (is.null(r1)) NA else r1$coefficients[, "exp(coef)"],
    if (is.null(r1)) NA else r1$coefficients[, "Pr(>|z|)"],
    if (is.null(r2)) NA else r2$coefficients["z", "exp(coef)"],
    if (is.null(r2)) NA else r2$coefficients["z", "Pr(>|z|)"])
}
m <- t(vapply(genes, one, numeric(4)))
colnames(m) <- c("HR_tcga", "p_tcga", "HR_adj", "p_adj")
res <- data.frame(gene = genes, m, row.names = NULL)
res <- merge(disc[, c("gene", "HR_disc", "FDR_disc", "logFC")], res, by = "gene")
res <- merge(res, comp, by = "gene")
res$sign <- ifelse(res$HR_disc > 1, 1, -1)
res$same_dir_uni <- (res$HR_tcga > 1) == (res$HR_disc > 1)
res$same_dir_adj <- (res$HR_adj > 1) == (res$HR_disc > 1)
res$replicated_uni <- !is.na(res$p_tcga) & res$p_tcga < 0.05 & res$same_dir_uni
res$replicated_adj <- !is.na(res$p_adj) & res$p_adj < 0.05 & res$same_dir_adj
res$delta_logHR <- res$sign * log(res$HR_tcga) - res$sign * log(res$HR_disc)   # effect change disc -> TCGA
res$class <- ifelse(res$share_bulk_stromal > 0.5, "stromal-dominant",
             ifelse(res$share_bulk_epithelial > 0.5, "epithelial-dominant",
             ifelse(res$share_bulk_immune > 0.5, "immune-dominant", "mixed")))
res <- res[!is.na(res$p_tcga), ]
write.csv(res, "systematic_replication.csv", row.names = FALSE)

## ---- summary and tests --------------------------------------------------------------------
sink("systematic_summary.txt", split = TRUE)
cat("Genes analysed (FDR_disc < 0.05, present in TCGA and in the single-cell table):", nrow(res), "\n")
cat("Overall replication, univariable:", sum(res$replicated_uni), "(",
    round(100 * mean(res$replicated_uni), 1), "% )\n")
cat("Overall replication after adjustment:", sum(res$replicated_adj), "(",
    round(100 * mean(res$replicated_adj), 1), "% )\n")
cat("Median risk-aligned log(HR): discovery", round(median(res$sign * log(res$HR_disc)), 3),
    "| TCGA", round(median(res$sign * log(res$HR_tcga)), 3), "\n\n")

tab <- do.call(rbind, lapply(split(res, res$class), function(x)
  data.frame(class = x$class[1], n = nrow(x), rep_uni = sum(x$replicated_uni),
             pct_uni = round(100 * mean(x$replicated_uni), 1), rep_adj = sum(x$replicated_adj),
             pct_adj = round(100 * mean(x$replicated_adj), 1),
             median_delta_logHR = round(median(x$delta_logHR), 3))))
print(tab, row.names = FALSE)

cat("\n--- PRE-SPECIFIED TEST: stromal-dominant vs epithelial-dominant, replication after adjustment ---\n")
two <- res[res$class %in% c("stromal-dominant", "epithelial-dominant"), ]
if (length(unique(two$class)) == 2) {
  tt <- table(two$class, two$replicated_adj)
  print(tt)
  if (all(dim(tt) == c(2, 2))) tryCatch(print(fisher.test(tt)), error = function(e) cat("test failed:", conditionMessage(e), "\n"))
  two$stromal <- as.integer(two$class == "stromal-dominant")
  two$abs_logHR_disc <- abs(log(two$HR_disc))
  mod <- tryCatch(glm(replicated_adj ~ stromal + abs_logHR_disc, data = two, family = binomial),
                  error = function(e) NULL)
  if (!is.null(mod)) { cat("\nLogistic regression (adjusted for discovery effect size):\n"); print(summary(mod)$coefficients) }
} else cat("Not enough genes in both classes for the test.\n")

cat("\n--- univariable replication, same comparison ---\n")
if (length(unique(two$class)) == 2) {
  tt2 <- table(two$class, two$replicated_uni); print(tt2)
  if (all(dim(tt2) == c(2, 2))) tryCatch(print(fisher.test(tt2)), error = function(e) cat("test failed:", conditionMessage(e), "\n"))
}
cat("\n--- continuous: stromal share vs change in effect size (discovery -> TCGA) ---\n")
tryCatch(print(cor.test(res$share_bulk_stromal, res$delta_logHR, method = "spearman", exact = FALSE)),
         error = function(e) cat("test failed:", conditionMessage(e), "\n"))
cat("\n--- continuous: stromal share vs replication after adjustment (Wilcoxon on share) ---\n")
tryCatch(print(wilcox.test(share_bulk_stromal ~ replicated_adj, data = res)),
         error = function(e) cat("test failed:", conditionMessage(e), "\n"))
sink()

pdf("systematic_replication.pdf", width = 7, height = 5.5)
tt3 <- tab[tab$n >= 5, ]
barplot(rbind(tt3$pct_uni, tt3$pct_adj), beside = TRUE, names.arg = paste0(tt3$class, "\n(n=", tt3$n, ")"),
        col = c("grey70", "firebrick"), ylab = "% of genes replicated in TCGA", cex.names = 0.8,
        main = "Replication by cellular compartment")
legend("topright", c("univariable", "adjusted (age, gender, stage)"), fill = c("grey70", "firebrick"), bty = "n")
dev.off()
cat("Done. See systematic_summary.txt, systematic_replication.csv, systematic_replication.pdf\n")
