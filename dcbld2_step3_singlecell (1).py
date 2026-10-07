"""
DCBLD2 expression by cell type in GSE132465 (single-cell CRC, 10X).
Reads the big count matrix in chunks (low memory), computes library size per
cell, then log-normalised DCBLD2 expression per cell type / subtype.

Put these two files (from the GSE132465 page, 'Supplementary file' section)
in the working folder and edit the two names below if yours differ:
    MATRIX  : raw UMI count matrix  (genes x cells, tab-separated, .gz)
    ANNOT   : cell annotation file  (cell id, patient, class, cell type, ...)

Run (Anaconda Prompt):   conda activate base
                         cd /d D:\\crc
                         python dcbld2_step3_singlecell.py

"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

MATRIX = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
ANNOT = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
GENE = "DCBLD2"
MARKERS = ["EPCAM", "KRT19", "CDX2", "COL1A1", "COL1A2", "DCN", "FAP", "PECAM1",
           "VWF", "PTPRC", "CD3E", "CD79A", "LYZ", "ST3GAL6", "TGFBR1", "TGFB1", "TGFB2"]
TARGETS = [GENE] + MARKERS
CHUNK = 300  # rows (genes) per chunk; lower it if you run out of memory

for f in (MATRIX, ANNOT):
    if not os.path.exists(f):
        sys.exit("Missing file: %s\nDownload it from the GSE132465 GEO page "
                 "(Supplementary file section) into this folder, or edit the "
                 "file name at the top of the script." % f)

# ---- annotation ----------------------------------------------------------
ann = pd.read_csv(ANNOT, sep="\t")
print("Annotation columns:", ann.columns.tolist())
print(ann.head(3).to_string())


def find_col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


idx_c = find_col(ann, "Index", "cell", "barcode") or ann.columns[0]
type_c = find_col(ann, "Cell_type", "celltype")
sub_c = find_col(ann, "Cell_subtype", "subtype")
cls_c = find_col(ann, "Class")
if type_c is None:
    sys.exit("Could not find a cell-type column. Send me the column list above.")
ann = ann.set_index(idx_c)

# ---- stream the count matrix ------------------------------------------------
lib = None
cells = None
rows = {}
reader = pd.read_csv(MATRIX, sep="\t", index_col=0, chunksize=CHUNK)
for i, ch in enumerate(reader):
    if lib is None:
        cells = ch.columns.to_numpy()
        lib = np.zeros(len(cells))
    lib += ch.to_numpy(dtype=np.float64).sum(axis=0)
    for g in TARGETS:
        if g in ch.index:
            r = ch.loc[g]
            if isinstance(r, pd.DataFrame):
                r = r.iloc[0]
            rows[g] = r.to_numpy(dtype=np.float64)
    if i % 10 == 0:
        print("read about %d genes ..." % ((i + 1) * CHUNK), flush=True)

if GENE not in rows:
    sys.exit("%s was not found in the matrix row names." % GENE)
print("Finished reading. Cells in matrix:", len(cells))

df = pd.DataFrame(index=cells)
df["lib"] = lib
for g, v in rows.items():
    df[g + "_count"] = v
    df[g + "_expr"] = np.log1p(v / np.maximum(lib, 1) * 1e4)
df = df.join(ann, how="inner")
print("Cells matched to annotation: %d of %d" % (len(df), len(cells)))
if len(df) == 0:
    sys.exit("No cell IDs matched between matrix and annotation. "
             "Send me a few cell IDs from each file.")


def summarize(by):
    g = df.groupby(by)
    out = pd.DataFrame({
        "n_cells": g.size(),
        "mean_logexpr": g[GENE + "_expr"].mean(),
        "pct_expressing": g[GENE + "_count"].apply(lambda s: (s > 0).mean() * 100),
    })
    return out.sort_values("mean_logexpr", ascending=False)


by_type = summarize(type_c)
by_type.to_csv("DCBLD2_sc_by_celltype.csv")
print("\n", by_type.round(3).to_string())

main = by_type
if sub_c:
    by_sub = summarize(sub_c)
    by_sub.to_csv("DCBLD2_sc_by_subtype.csv")
    print("\n", by_sub.round(3).to_string())
    main = by_sub
if cls_c:
    summarize([type_c, cls_c]).to_csv("DCBLD2_sc_by_celltype_and_class.csv")

# sanity check that the annotation matches known markers
mk = [m + "_expr" for m in MARKERS if m in rows]
df.groupby(type_c)[mk].mean().round(3).to_csv("DCBLD2_sc_marker_check.csv")

# ---- plot -----------------------------------------------------------------------
m = main.sort_values("mean_logexpr")
fig, ax = plt.subplots(1, 2, figsize=(11, max(4, 0.35 * len(m) + 1)), sharey=True)
ax[0].barh(m.index.astype(str), m["mean_logexpr"], color="firebrick")
ax[0].set_xlabel("Mean log-normalised %s" % GENE)
ax[1].barh(m.index.astype(str), m["pct_expressing"], color="steelblue")
ax[1].set_xlabel("% of cells expressing %s" % GENE)
fig.suptitle("%s in GSE132465 by %s" % (GENE, sub_c or type_c))
fig.tight_layout()
fig.savefig("DCBLD2_sc_barplot.png", dpi=200)
print("\nDone. Saved DCBLD2_sc_by_celltype.csv, DCBLD2_sc_marker_check.csv, "
      "DCBLD2_sc_barplot.png (and subtype/class tables if available).")
