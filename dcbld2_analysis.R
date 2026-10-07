# =====================================================================
# DCBLD2 characterization in colorectal cancer using GEO cohorts only
#   Discovery : GSE39582 (GPL570)   Validation: GSE17536 (GPL570)
# Needs in the working folder (already downloaded):
#   GSE39582_series_matrix.txt.gz   GSE17536_series_matrix.txt.gz
# Outputs:
#   DCBLD2_figures.pdf              (all plots, one per page)
#   DCBLD2_cox_results.csv          (univariable + multivariable Cox)
#   DCBLD2_marker_correlations.csv  (stromal/epithelial marker check)
#   DCBLD2_correlated_genes.csv     (co-expression, GSE39582 tumors)
#   DCBLD2_GO_positive.csv / DCBLD2_GO_negative.csv
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

GENE     <- "DCBLD2"
DISC_ACC <- "GSE39582"
VAL_ACC  <- "GSE17536"

## ---- 0. Packages ----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("GEOquery", "limma", "Biobase", "AnnotationDbi",
            "hgu133plus2.db", "org.Hs.eg.db", "GO.db"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")

suppressMessages({
  library(GEOquery); library(limma); library(Biobase)
  library(AnnotationDbi); library(hgu133plus2.db); library(org.Hs.eg.db)
  library(survival)
})

## ---- Helpers ---------------------------------------------------------
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

zscore_mean <- function(ex, genes, samples) {
  g <- intersect(genes, rownames(ex))
  m <- ex[g, samples, drop = FALSE]
  colMeans(t(scale(t(m))), na.rm = TRUE)
}

fib   <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "FAP", "PDGFRB", "LUM")
endo  <- c("PECAM1", "VWF", "CDH5", "KDR", "CLDN5")
epi   <- c("EPCAM", "KRT19", "CDH1", "VIL1", "CDX2")

cox_table <- function(m, label, model) {
  s <- summary(m)
  data.frame(cohort = label, model = model, term = rownames(s$coefficients),
             HR = s$coefficients[, "exp(coef)"],
             lower95 = s$conf.int[, "lower .95"], upper95 = s$conf.int[, "upper .95"],
             p = s$coefficients[, "Pr(>|z|)"], n = m$n, events = m$nevent,
             row.names = NULL)
}

run_cox <- function(ex, os, pd, label, stromal = NULL) {
  d <- os[os$sample %in% colnames(ex), ]
  d$gene <- as.numeric(scale(ex[GENE, d$sample]))
  out <- list(cox_table(coxph(Surv(time, event) ~ gene, data = d), label, "univariable"))
  cv <- make_cov(pd, d$sample, label)
  if (ncol(cv) > 1) d <- cbind(d, cv[, -1, drop = FALSE])
  covs <- setdiff(colnames(cv), "sample")
  covs <- covs[vapply(covs, function(k)
    mean(is.na(d[[k]])) < 0.3 && length(unique(na.omit(d[[k]]))) > 1, logical(1))]
  if (length(covs)) {
    f <- as.formula(paste("Surv(time, event) ~ gene +", paste(covs, collapse = " + ")))
    out[[2]] <- cox_table(coxph(f, data = d), label, paste("multivariable:", paste(covs, collapse = "+")))
  }
  if (!is.null(stromal)) {
    d$stromal <- as.numeric(scale(stromal[d$sample]))
    out[[length(out) + 1]] <- cox_table(coxph(Surv(time, event) ~ gene + stromal, data = d),
                                        label, "adjusted for stromal score")
  }
  do.call(rbind, out)
}

km_plot <- function(ex, os, title) {
  d <- os[os$sample %in% colnames(ex), ]
  d$grp <- factor(ifelse(ex[GENE, d$sample] > median(ex[GENE, d$sample]), "High", "Low"),
                  levels = c("Low", "High"))
  fit <- survfit(Surv(time, event) ~ grp, data = d)
  lr <- survdiff(Surv(time, event) ~ grp, data = d)
  p <- 1 - pchisq(lr$chisq, 1)
  plot(fit, col = c("steelblue", "firebrick"), lwd = 2, mark.time = TRUE,
       xlab = "Time (months)", ylab = "Overall survival", main = title)
  legend("bottomleft", c(paste(GENE, "low"), paste(GENE, "high")),
         col = c("steelblue", "firebrick"), lwd = 2, bty = "n")
  legend("topright", sprintf("log-rank p = %.2g (median split)", p), bty = "n")
}

## ---- 1. Load both cohorts ------------------------------------------------
g1 <- load_gse(DISC_ACC); pd1 <- pData(g1); ex1 <- to_gene(g1)
g2 <- load_gse(VAL_ACC);  pd2 <- pData(g2); ex2 <- to_gene(g2)
stopifnot(GENE %in% rownames(ex1), GENE %in% rownames(ex2))

isnorm <- apply(pd1, 1, function(r) any(grepl("non.?tumou?r", r, ignore.case = TRUE)))
tum1 <- colnames(ex1)[!isnorm]
cat("GSE39582: normal =", sum(isnorm), " tumor =", length(tum1), "\n")

os1 <- get_os(pd1, DISC_ACC); os1 <- os1[os1$sample %in% tum1, ]
os2 <- get_os(pd2, VAL_ACC)

## ---- 2. Cox models (univariable, multivariable, stromal-adjusted) ---------
st1 <- zscore_mean(ex1, c(fib, endo), tum1)
st2 <- zscore_mean(ex2, c(fib, endo), colnames(ex2))
cox_res <- rbind(run_cox(ex1, os1, pd1, DISC_ACC, stromal = st1),
                 run_cox(ex2, os2, pd2, VAL_ACC,  stromal = st2))
write.csv(cox_res, "DCBLD2_cox_results.csv", row.names = FALSE)
print(cox_res[cox_res$term == "gene", c("cohort", "model", "HR", "lower95", "upper95", "p")])

## ---- 3. Marker correlations (is DCBLD2 a stromal signal?) ------------------
mk <- data.frame(gene = c(fib, endo, epi),
                 compartment = c(rep("Fibroblast", length(fib)),
                                 rep("Endothelial", length(endo)),
                                 rep("Epithelial", length(epi))))
mk <- mk[mk$gene %in% rownames(ex1), ]
x <- ex1[GENE, tum1]
ct <- lapply(mk$gene, function(g) {
  r <- suppressWarnings(cor.test(x, ex1[g, tum1], method = "spearman"))
  c(rho = unname(r$estimate), p = r$p.value)
})
mk <- cbind(mk, do.call(rbind, ct))
write.csv(mk, "DCBLD2_marker_correlations.csv", row.names = FALSE)

## ---- 4. Figures -------------------------------------------------------------
pdf("DCBLD2_figures.pdf", width = 7, height = 6)

# 4a tumor vs normal (GSE39582)
grp <- factor(ifelse(isnorm, "Normal", "Tumor"), levels = c("Normal", "Tumor"))
pw <- wilcox.test(ex1[GENE, ] ~ grp)$p.value
boxplot(ex1[GENE, ] ~ grp, col = c("grey80", "salmon"), ylab = paste(GENE, "log2 expression"),
        main = paste0(GENE, " in ", DISC_ACC), outline = FALSE)
stripchart(ex1[GENE, ] ~ grp, vertical = TRUE, method = "jitter", pch = 16,
           col = rgb(0, 0, 0, 0.3), add = TRUE)
legend("topleft", sprintf("Wilcoxon p = %.2g", pw), bty = "n")

# 4b Kaplan-Meier in both cohorts
km_plot(ex1[, tum1], os1, paste(GENE, "-", DISC_ACC))
km_plot(ex2, os2, paste(GENE, "-", VAL_ACC))

# 4c stage association (discovery tumors)
cv1 <- make_cov(pd1, tum1, DISC_ACC)
if ("stage" %in% colnames(cv1) && sum(!is.na(cv1$stage)) > 30) {
  v <- ex1[GENE, tum1]
  p <- wilcox.test(v ~ cv1$stage)$p.value
  boxplot(v ~ cv1$stage, col = c("grey85", "orange"), outline = FALSE,
          ylab = paste(GENE, "log2 expression"), main = paste(GENE, "by stage -", DISC_ACC))
  legend("topleft", sprintf("Wilcoxon p = %.2g", p), bty = "n")
}

# 4d marker correlations
mk2 <- mk[order(mk$rho), ]
cols <- c(Fibroblast = "darkorange", Endothelial = "purple", Epithelial = "darkgreen")[mk2$compartment]
par(mar = c(4, 6, 3, 1))
barplot(mk2$rho, names.arg = mk2$gene, horiz = TRUE, las = 1, col = cols,
        xlab = paste("Spearman rho with", GENE, "(tumors,", DISC_ACC, ")"),
        main = "Stromal vs epithelial marker correlation")
legend("bottomright", names(table(mk2$compartment)),
       fill = c(Fibroblast = "darkorange", Endothelial = "purple", Epithelial = "darkgreen")[names(table(mk2$compartment))],
       bty = "n")
dev.off()
cat("Saved DCBLD2_figures.pdf, DCBLD2_cox_results.csv, DCBLD2_marker_correlations.csv\n")

## ---- 5. Co-expression + GO enrichment (limma::goana) -----------------------
tryCatch({
  cm <- cor(t(ex1[rownames(ex1) != GENE, tum1]), ex1[GENE, tum1], method = "spearman")
  cdf <- data.frame(gene = rownames(cm), rho = cm[, 1], row.names = NULL)
  cdf <- cdf[order(-cdf$rho), ]
  write.csv(cdf, "DCBLD2_correlated_genes.csv", row.names = FALSE)

  to_entrez <- function(s) unique(na.omit(unname(
    suppressMessages(mapIds(org.Hs.eg.db, keys = s, column = "ENTREZID",
                            keytype = "SYMBOL", multiVals = "first")))))
  univ <- to_entrez(cdf$gene)
  for (side in c("positive", "negative")) {
    top <- if (side == "positive") head(cdf$gene, 300) else tail(cdf$gene, 300)
    res <- goana(to_entrez(top), universe = univ, species = "Hs")
    res$FDR <- p.adjust(res$P.DE, "BH")
    out <- topGO(res, ontology = "BP", number = 50)
    write.csv(out, sprintf("DCBLD2_GO_%s.csv", side))
    cat("Saved DCBLD2_GO_", side, ".csv\n", sep = "")
  }
}, error = function(e) message("Co-expression/GO step failed: ", conditionMessage(e)))

cat("Done.\n")
