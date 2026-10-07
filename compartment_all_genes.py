"""
Cellular compartment of origin for EVERY gene in the GSE132465 single-cell CRC data.

One pass over the count matrix. Only TUMOR-class cells are used. For each gene:
  share_bulk_<compartment>  : fraction of the gene's UMIs coming from that compartment
                              (abundance-weighted, i.e. what a bulk RNA-seq sample "sees")
  share_spec_<compartment>  : cell-type-normalised share (CPM per compartment, then
                              normalised) = expression specificity independent of abundance
Compartments: Epithelial, Stromal, Immune (T, B, myeloid, mast cells).

Run (Anaconda Prompt):  conda activate base
                        cd /d D:\\crc
                        python compartment_all_genes.py > comp_log.txt 2>&1
Needs the two GSE132465 files already in the folder.
Output: compartment_all_genes.csv  (~10-20 min; the window stays silent)
"""
import os
import sys
import numpy as np
import pandas as pd

MATRIX = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
ANNOT = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
CLASS_KEEP = "Tumor"
MIN_UMI = 30        # genes with fewer UMIs in tumor cells are dropped (too noisy)
CHUNK = 300

for f in (MATRIX, ANNOT):
    if not os.path.exists(f):
        sys.exit("Missing file: " + f)


def find_col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


ann = pd.read_csv(ANNOT, sep="\t")
print("Annotation columns:", ann.columns.tolist(), flush=True)
idx_c = find_col(ann, "Index", "cell", "barcode") or ann.columns[0]
type_c = find_col(ann, "Cell_type", "celltype")
cls_c = find_col(ann, "Class")
if type_c is None or cls_c is None:
    sys.exit("Need Cell_type and Class columns. Send me the column list above.")
ann = ann.set_index(idx_c)
tumor = ann[ann[cls_c] == CLASS_KEEP]
types = sorted(tumor[type_c].astype(str).unique())
print("Tumor cells: %d | cell types: %s" % (len(tumor), types), flush=True)


def compartment(t):
    t = t.lower()
    if "epithelial" in t:
        return "epithelial"
    if "stromal" in t:
        return "stromal"
    return "immune"


comp_of_type = [compartment(t) for t in types]
print("Compartment mapping:", dict(zip(types, comp_of_type)), flush=True)

lib_all, ind, sums, names = None, None, [], []
for i, ch in enumerate(pd.read_csv(MATRIX, sep="\t", index_col=0, chunksize=CHUNK)):
    if ind is None:
        cells = pd.Index(ch.columns)
        lib_all = np.zeros(len(cells))
        ind = np.zeros((len(cells), len(types)))
        pos = cells.get_indexer(tumor.index)
        code = pd.Categorical(tumor[type_c].astype(str), categories=types).codes
        ok = pos >= 0
        ind[pos[ok], code[ok]] = 1.0
        print("Tumor cells found in matrix: %d" % ok.sum(), flush=True)
        if ok.sum() == 0:
            sys.exit("No annotated tumor cells match the matrix cell IDs.")
    vals = ch.to_numpy(dtype=np.float64)
    lib_all += vals.sum(axis=0)
    sums.append(vals @ ind)
    names.extend(ch.index.tolist())
    if i % 10 == 0:
        print("read about %d genes ..." % ((i + 1) * CHUNK), flush=True)

counts = pd.DataFrame(np.vstack(sums), index=names, columns=types)
counts = counts[~counts.index.duplicated(keep="first")]
lib_type = pd.Series(lib_all @ ind, index=types)
total = counts.sum(axis=1)
counts = counts[total >= MIN_UMI]
total = counts.sum(axis=1)
cpm = counts.div(lib_type, axis=1) * 1e6

out = pd.DataFrame({"gene": counts.index, "total_tumor_UMI": total.values})
for c in ("epithelial", "stromal", "immune"):
    cols = [t for t, k in zip(types, comp_of_type) if k == c]
    out["share_bulk_" + c] = (counts[cols].sum(axis=1) / total).values
    spec_num = cpm[cols].sum(axis=1)
    out["share_spec_" + c] = (spec_num / cpm.sum(axis=1)).values
for t in types:
    out["UMI_" + t.replace(" ", "_")] = counts[t].values
out["dominant_bulk"] = out[["share_bulk_epithelial", "share_bulk_stromal", "share_bulk_immune"]] \
    .idxmax(axis=1).str.replace("share_bulk_", "", regex=False)
out.to_csv("compartment_all_genes.csv", index=False)
print("\nGenes kept: %d" % len(out))
print("Dominant compartment (bulk share):")
print(out["dominant_bulk"].value_counts().to_string())
print("Genes with >50%% of UMIs in one compartment: epithelial %d | stromal %d | immune %d" % (
    (out.share_bulk_epithelial > .5).sum(), (out.share_bulk_stromal > .5).sum(),
    (out.share_bulk_immune > .5).sum()))
for g in ("DCBLD2", "TGFB2", "CERCAM", "ISM1", "KRT17", "LAMC2", "KLK6", "TGFBRAP1"):
    r = out[out.gene == g]
    if len(r):
        r = r.iloc[0]
        print("%-9s epithelial %.2f | stromal %.2f | immune %.2f  (UMI %d)" % (
            g, r.share_bulk_epithelial, r.share_bulk_stromal, r.share_bulk_immune, r.total_tumor_UMI))
print("\nDone. Saved compartment_all_genes.csv")
