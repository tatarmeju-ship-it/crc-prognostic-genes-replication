# Do prognostic genes from bulk colorectal cancer transcriptomes replicate?
A test of cellular origin using single-cell data

Code and gene-level results for the manuscript of the same title ([journal / DOI to be added]).

## Summary

We identified genes associated with overall survival in a discovery cohort of colon cancer
(GSE39582) and tested replication in two independent cohorts (TCGA-COADREAD and GSE17536),
with univariable and age-, sex- and stage-adjusted Cox models. The cellular origin of each
gene's tumour expression was estimated from a public single-cell dataset (GSE132465).
A hypothesis that stromal-dominant genes replicate more often was generated in an exploratory
set of 223 genes and tested in a separate confirmatory set of 78 genes. DCBLD2 is used as a
worked example (single-cell expression, survival, Mendelian randomisation).

Key results:
- In pooled replication data, 31% and 43% of the discovery effect was retained (two gene sets).
- Only 8% and 18% of genes remained significant after adjustment for age, sex and stage.
- The primary hypothesis was not confirmed at the 0.05 level (Fisher p = 0.062).

## Status of the analyses

| Analysis | Status |
|---|---|
| Exploratory set (223 genes) | Hypothesis-generating |
| Confirmatory set (78 genes), primary test | Planned before execution (see below) |
| Sensitivity and post hoc analyses | Secondary / exploratory, uncorrected |
| DCBLD2 worked example, Mendelian randomisation | Descriptive |

**Analysis plan.** The hypothesis, gene-selection and exclusion rules, primary outcome and
primary test of the confirmatory analysis are written in the header of
`scripts/04_confirmatory/confirmatory_analysis.R`. The plan was not formally registered and
this repository was created after the analyses were run, so commit dates do not verify when the
plan was written. We therefore do not claim pre-registration.

## Repository contents

| Path | Description |
|---|---|
| `scripts/` | Analysis scripts, in the order of the table below |
| `results/` | Gene-level results (confirmatory 78 genes, exploratory 223 genes) |
| `supplementary_tables/` | S1-S3 Tables (Excel) |
| `LICENSE` | MIT (code) |

## Run order

| Step | Script | Output |
|---|---|---|
| 1 | `crc_geo_candidates.R` | `candidates_discovery.csv` |
| 2 | `compartment_all_genes.py` | `compartment_all_genes.csv` |
| 3 | `systematic_replication.R`, `systematic_replication2.R` | exploratory results |
| 4 | `confirmatory_analysis.R`, `flow_numbers.R` | `confirmatory_results.csv`, gene-flow numbers |
| 5 | DCBLD2 scripts: `dcbld2_analysis.R`, `dcbld2_tcga_validation.R`, `gene_pseudobulk.py`, `extract_dcbld2_eqtl.py`, `extract_gwas_region.py`, `dcbld2_step4_mr.py` | worked-example results |
| 6 | `make_figures.py` | manuscript figures |

Scripts read their inputs from the working folder (file names are given in each script header).
They were written and run on Windows with Anaconda.

## Data sources (publicly available; not included here)

- Discovery: GSE39582 (GEO)
- Replication: TCGA-COADREAD (UCSC Xena), GSE17536 (GEO)
- Single-cell: GSE132465 (GEO)
- Genetic instruments: eQTLGen cis-eQTL summary statistics
- GWAS outcome: GCST90255675 (academic non-commercial licence; not redistributed)

All data are de-identified; no new patient data were generated.

## Software

R [version] with survival, limma, GEOquery, AnnotationDbi, hgu133plus2.db.
Python [version] with pandas, SciPy, matplotlib. Versions: `sessionInfo.txt`, `environment.yml`.

## AI assistance

Parts of the code and the manuscript were developed with the help of an AI language model
(Claude, Anthropic). All analyses were run and all results were reviewed by the authors.

## Citation

[Authors]. [Title]. [Journal, year]. Code archive: https://doi.org/[Zenodo DOI]
(see `CITATION.cff`).
