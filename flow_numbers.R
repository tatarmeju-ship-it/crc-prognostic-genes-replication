# =====================================================================
# Gene-flow numbers for the Methods / flow diagram (no new hypothesis testing).
# Repeats only the discovery step of confirmatory_analysis.R (a few minutes) and
# counts how many genes remain after each filter.  Output: flow_numbers.txt
# Needs: candidates_discovery.csv, compartment_all_genes.csv, GSE39582 and GSE17536
#        series matrix files, TCGA HiSeqV2 file.
# Run:   Rscript flow_numbers.R > flow_log.txt 2>&1
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

comp <- read.csv("compartment_all_genes.csv", stringsAsFactors = FALSE)
old_set <- read.csv("candidates_discovery.csv", stringsAsFactors = FALSE)
excluded <- old_set$gene[!is.na(old_set$FDR_disc) & old_set$FDR_disc < 0.05]
message("Genes excluded (already examined in the exploratory analysis): ", length(excluded))

## ---- discovery: univariable Cox for ALL genes in GSE39582 tumors ----------------------------
g1d <- load_gse("GSE39582"); pd1 <- pData(g1d); ex1 <- to_gene(g1d)
isnorm <- apply(pd1, 1, function(r) any(grepl("non.?tumou?r", r, ignore.case = TRUE)))
os1 <- get_os(pd1, "GSE39582"); os1 <- os1[os1$sample %in% colnames(ex1)[!isnorm], ]
message("GSE39582 tumors with OS: ", nrow(os1), ", deaths: ", sum(os1$event),
        " | fitting Cox models for ", nrow(ex1), " genes (a few minutes) ...")
md <- t(vapply(rownames(ex1), function(g) cox_est(ex1[g, os1$sample], os1, character(0))[1:3], numeric(3)))
colnames(md) <- c("b_disc", "se_disc", "p_disc")
disc <- data.frame(gene = rownames(ex1), md, row.names = NULL)
disc <- disc[!is.na(disc$p_disc), ]
disc$FDR_disc <- p.adjust(disc$p_disc, "BH")
disc$HR_disc <- exp(disc$b_disc)
message("Genes with FDR < 0.05 in discovery: ", sum(disc$FDR_disc < 0.05))

## ---- count genes after each filter ---------------------------------------------------------------
n_tested <- nrow(disc)
sig <- disc[disc$FDR_disc < 0.05, ]
s2  <- sig[!(sig$gene %in% excluded), ]
s3  <- s2[s2$gene %in% comp$gene, ]
gx  <- read.delim(gzfile(find1("HiSeqV2")), check.names = FALSE, stringsAsFactors = FALSE)[[1]]
s4  <- s3[s3$gene %in% gx, ]
g17 <- rownames(to_gene(load_gse("GSE17536")))
s5  <- s4[s4$gene %in% g17, ]
txt <- c(
  paste("Genes with a valid Cox model in the discovery cohort (GSE39582 tumors):", n_tested),
  paste("Genes with FDR < 0.05 in discovery:", nrow(sig)),
  paste("  of which already examined in the exploratory analysis (excluded):", sum(sig$gene %in% excluded)),
  paste("New genes after exclusion:", nrow(s2)),
  paste("  present in the single-cell compartment table:", nrow(s3)),
  paste("  also present in TCGA-COADREAD:", nrow(s4)),
  paste("  also present in GSE17536 (final confirmatory set, before dropping genes with failed models):", nrow(s5)),
  "The analysed set reported earlier was 78 genes.")
writeLines(txt, "flow_numbers.txt")
cat(txt, sep = "\n")
