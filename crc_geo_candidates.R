# =====================================================================
# CRC candidate genes using GEO only (no GDC/TCGA needed)
#   Discovery : GSE39582 (GPL570) -> tumor vs normal (limma) -> Cox (OS)
#   Validation: VAL_ACC (default GSE17538, GPL570) -> Cox (OS)
# Output: candidates_discovery.csv, candidates_ranked.csv
# =====================================================================
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)

DISC_ACC <- "GSE39582"
VAL_ACC  <- "GSE17536"   # must be a GPL570 cohort with overall-survival data

## ---- 0. Packages ----------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (p in c("GEOquery", "limma", "Biobase", "AnnotationDbi", "hgu133plus2.db"))
  if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
for (p in c("survival", "rentrez"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

suppressMessages({
  library(GEOquery); library(limma); library(Biobase)
  library(AnnotationDbi); library(hgu133plus2.db)
  library(survival); library(rentrez)
})

## ---- Helpers ---------------------------------------------------------
# Uses a local <ACC>_series_matrix.txt.gz if present (manual download),
# otherwise downloads from GEO.
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

# probe -> gene symbol (GPL570), one probe per gene (highest mean)
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

# find overall-survival columns automatically and report which ones are used
get_os <- function(pd, label) {
  cn <- colnames(pd)
  tcol <- grep("(^|[^a-z])os.*(delay|time|month)|overall.*(time|month|delay|surviv)", cn,
               ignore.case = TRUE, value = TRUE)
  ecol <- grep("(^|[^a-z])os.*(event|status)|overall.*(event|status)|vital|death", cn,
               ignore.case = TRUE, value = TRUE)
  if (!length(tcol) || !length(ecol)) {
    message("Columns in ", label, ":\n", paste(cn, collapse = "\n"))
    stop("Could not find OS time/event columns in ", label,
         ". Send me the column list printed above.")
  }
  message(label, ": OS time column  = ", tcol[1])
  message(label, ": OS event column = ", ecol[1])
  message(label, ": example values  = ", paste(head(unique(as.character(pd[[ecol[1]]])), 5), collapse = " | "))
  tm <- suppressWarnings(as.numeric(as.character(pd[[tcol[1]]])))
  ev <- suppressWarnings(as.integer(as.character(pd[[ecol[1]]])))
  d <- data.frame(sample = rownames(pd), time = tm, event = ev)
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
  if (nrow(d) < 50) stop("Fewer than 50 patients with usable OS in ", label,
                         ". Check the printed columns/values and send them to me.")
  message(label, ": patients with OS = ", nrow(d), ", events = ", sum(d$event))
  d
}

cox_one <- function(x, d) {
  d$x <- as.numeric(scale(x))
  s <- tryCatch(summary(coxph(Surv(time, event) ~ x, data = d)), error = function(e) NULL)
  if (is.null(s)) return(c(NA_real_, NA_real_))
  c(s$coefficients[, "exp(coef)"], s$coefficients[, "Pr(>|z|)"])
}

run_cox <- function(ex, genes, d, suffix) {
  genes <- intersect(genes, rownames(ex))
  m <- t(apply(ex[genes, d$sample, drop = FALSE], 1, cox_one, d = d))
  out <- data.frame(gene = genes, m, row.names = NULL)
  colnames(out) <- c("gene", paste0("HR_", suffix), paste0("p_", suffix))
  out
}
get_os <- function(pd, label) {
  cn <- colnames(pd)
  tcol <- grep("overall.*(time|month|delay)|(^|[^a-z])os.*(delay|time|month)",
               cn, ignore.case = TRUE, value = TRUE)
  ecol <- grep("overall.*event|(^|[^a-z])os.*event|any cause",
               cn, ignore.case = TRUE, value = TRUE)
  used_dss <- FALSE
  if (!length(ecol)) {
    ecol <- grep("event|status|vital|death", cn, ignore.case = TRUE, value = TRUE)
    used_dss <- TRUE
  }
  if (!length(tcol) || !length(ecol)) {
    message("Columns in ", label, ":\n", paste(cn, collapse = "\n"))
    stop("Could not find survival columns in ", label, ". Send me the list above.")
  }
  message(label, ": time column  = ", tcol[1])
  message(label, ": event column = ", ecol[1])
  if (used_dss) message(label, ": WARNING - no overall-survival event column; using this one instead")
  raw <- tolower(trimws(as.character(pd[[ecol[1]]])))
  message(label, ": event values = ", paste(head(unique(raw), 6), collapse = " | "))
  ev <- ifelse(grepl("^(no death|alive|no|0|censored|false)", raw), 0L,
        ifelse(grepl("^(death|dead|deceased|yes|1|true)", raw), 1L, NA_integer_))
  tm <- suppressWarnings(as.numeric(as.character(pd[[tcol[1]]])))
  d <- data.frame(sample = rownames(pd), time = tm, event = ev)
  d <- d[!is.na(d$time) & d$time > 0 & !is.na(d$event), ]
  if (nrow(d) < 50) stop("Fewer than 50 patients with usable survival in ", label,
                         ". Send me the messages printed above.")
  message(label, ": patients = ", nrow(d), ", events = ", sum(d$event))
  d
}
## ---- 1. Discovery cohort: tumor vs normal ----------------------------
g1  <- load_gse(DISC_ACC)
pd1 <- pData(g1)
ex1 <- to_gene(g1)
cat("Genes in discovery matrix:", nrow(ex1), "\n")

isnorm <- apply(pd1, 1, function(r) any(grepl("non.?tumou?r", r, ignore.case = TRUE)))
cat("Normal samples:", sum(isnorm), " Tumor samples:", sum(!isnorm), "\n")
if (sum(isnorm) < 5) {
  print(head(pd1[, grep("tissue|title|source", colnames(pd1), ignore.case = TRUE), drop = FALSE], 20))
  stop("Could not identify normal samples. Send me the table printed above.")
}
grp <- factor(ifelse(isnorm, "Normal", "Tumor"), levels = c("Normal", "Tumor"))
fit <- eBayes(lmFit(ex1, model.matrix(~ grp)))
tt  <- topTable(fit, coef = 2, number = Inf)
tt$gene <- rownames(tt)
deg <- subset(tt, adj.P.Val < 0.01 & logFC > 0.585)   # up in tumor, >1.5-fold
cat("Up-regulated DEGs:", nrow(deg), "\n")

## ---- 2. Discovery cohort: Cox in tumors -------------------------------
os1 <- get_os(pd1, DISC_ACC)
os1 <- os1[os1$sample %in% colnames(ex1)[!isnorm], ]
cx1 <- run_cox(ex1, deg$gene, os1, "disc")
cx1$FDR_disc <- p.adjust(cx1$p_disc, "BH")
disc <- merge(cx1, deg[, c("gene", "logFC", "adj.P.Val")], by = "gene")
disc <- disc[order(disc$FDR_disc), ]
write.csv(disc, "candidates_discovery.csv", row.names = FALSE)
cat("Saved candidates_discovery.csv (", nrow(disc), " genes)\n")

## ---- 3. Validation cohort ----------------------------------------------
g2  <- load_gse(VAL_ACC)
pd2 <- pData(g2)
ex2 <- to_gene(g2)
os2 <- get_os(pd2, VAL_ACC)
os2 <- os2[os2$sample %in% colnames(ex2), ]
cx2 <- run_cox(ex2, disc$gene, os2, "val")

## ---- 4. Merge, filter, rank ---------------------------------------------
all <- merge(disc, cx2, by = "gene")
all$same_direction <- (all$HR_disc > 1) == (all$HR_val > 1)
final <- subset(all, FDR_disc < 0.05 & p_val < 0.05 & same_direction)
final <- final[order(final$FDR_disc + final$p_val), ]
cat("Candidates passing both cohorts:", nrow(final), "\n")

## ---- 5. PubMed novelty check (top 100) ------------------------------------
top <- head(final, 100)
if (nrow(top) > 0) {
  top$pubmed_CRC_hits <- vapply(top$gene, function(g) {
    Sys.sleep(0.4)
    tryCatch(entrez_search(db = "pubmed", retmax = 0,
      term = sprintf('%s[Title/Abstract] AND (colorectal OR colon OR rectal)', g))$count,
      error = function(e) NA_integer_)
  }, integer(1))
}
write.csv(top, "candidates_ranked.csv", row.names = FALSE)
cat("Done. Open candidates_ranked.csv\n")
