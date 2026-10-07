"""
Patient-level (pseudobulk) check of a gene in the GSE132465 single-cell CRC data.

For each cell type, counts are summed per patient and tissue class (Tumor/Normal),
converted to log2 CPM, and tumor vs normal is tested across PATIENTS (paired
Wilcoxon where a patient has both), which avoids treating thousands of cells
from the same person as independent.

Run (Anaconda Prompt):   conda activate base
                         cd /d D:\\crc
                         python gene_pseudobulk.py DCBLD2
Needs the same two GSE132465 files as gene_singlecell.py in the folder.
Outputs: <GENE>_sc_pseudobulk_by_patient.csv, <GENE>_sc_pseudobulk_tests.csv,
         <GENE>_sc_pseudobulk.png
NOTE: reading the matrix takes about as long as before (~10-20 min).
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENE = sys.argv[1] if len(sys.argv) > 1 else "DCBLD2"
MATRIX = "GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz"
ANNOT = "GSE132465_GEO_processed_CRC_10X_cell_annotation.txt.gz"
MIN_CELLS = 20      # minimum cells per patient x class x cell type
CHUNK = 300

for f in (MATRIX, ANNOT):
    if not os.path.exists(f):
        sys.exit("Missing file: %s (same files as for gene_singlecell.py)" % f)


def find_col(df, *names):
    low = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


ann = pd.read_csv(ANNOT, sep="\t")
print("Annotation columns:", ann.columns.tolist())
idx_c = find_col(ann, "Index", "cell", "barcode") or ann.columns[0]
type_c = find_col(ann, "Cell_type", "celltype")
cls_c = find_col(ann, "Class")
pat_c = find_col(ann, "Patient", "patient_id", "Subject")
samp_c = find_col(ann, "Sample")
if type_c is None or cls_c is None:
    sys.exit("Need Cell_type and Class columns. Send me the column list printed above.")
ann = ann.set_index(idx_c)
if pat_c is None:
    if samp_c is None:
        sys.exit("No Patient or Sample column found. Send me the column list printed above.")
    ann["Patient"] = ann[samp_c].astype(str).str.extract(r"^([A-Za-z]+\d+)")[0]
    pat_c = "Patient"
    print("Patient ID derived from the Sample column.")

# ---- stream the matrix: only library size and the target gene are kept ----
lib, cells, gene_counts = None, None, None
for i, ch in enumerate(pd.read_csv(MATRIX, sep="\t", index_col=0, chunksize=CHUNK)):
    if lib is None:
        cells = ch.columns.to_numpy()
        lib = np.zeros(len(cells))
    lib += ch.to_numpy(dtype=np.float64).sum(axis=0)
    if GENE in ch.index:
        r = ch.loc[GENE]
        if isinstance(r, pd.DataFrame):
            r = r.iloc[0]
        gene_counts = r.to_numpy(dtype=np.float64)
    if i % 10 == 0:
        print("read about %d genes ..." % ((i + 1) * CHUNK), flush=True)
if gene_counts is None:
    sys.exit("%s was not found in the matrix." % GENE)

df = pd.DataFrame({"lib": lib, "count": gene_counts}, index=cells).join(ann, how="inner")
df = df[df[cls_c].isin(["Tumor", "Normal"])]
print("Cells used (Tumor/Normal only): %d" % len(df))
df["expr_cell"] = (df["count"] > 0).astype(float)

g = (df.groupby([type_c, pat_c, cls_c])
       .agg(n_cells=("lib", "size"), counts=("count", "sum"), lib=("lib", "sum"),
            frac_expressing=("expr_cell", "mean"))
       .reset_index())
g["log2cpm"] = np.log2(g["counts"] / g["lib"] * 1e6 + 1)
g = g[g["n_cells"] >= MIN_CELLS]
g.to_csv("%s_sc_pseudobulk_by_patient.csv" % GENE, index=False)

rows = []
for t, sub in g.groupby(type_c):
    piv = sub.pivot(index=pat_c, columns=cls_c, values="log2cpm")
    if not {"Tumor", "Normal"} <= set(piv.columns):
        continue
    both = piv.dropna(subset=["Tumor", "Normal"])
    allt = piv["Tumor"].dropna()
    alln = piv["Normal"].dropna()
    row = dict(cell_type=t, patients_tumor=len(allt), patients_normal=len(alln), paired_patients=len(both),
               mean_log2cpm_tumor=allt.mean(), mean_log2cpm_normal=alln.mean())
    if len(both) >= 4 and (both["Tumor"] - both["Normal"]).abs().sum() > 0:
        row["paired_wilcoxon_p"] = stats.wilcoxon(both["Tumor"], both["Normal"]).pvalue
        row["mean_paired_diff"] = (both["Tumor"] - both["Normal"]).mean()
    if len(allt) >= 3 and len(alln) >= 3:
        row["unpaired_mannwhitney_p"] = stats.mannwhitneyu(allt, alln).pvalue
    rows.append(row)
res = pd.DataFrame(rows).sort_values("paired_patients", ascending=False)
res.to_csv("%s_sc_pseudobulk_tests.csv" % GENE, index=False)
pd.set_option("display.width", 220)
print(res.round(4).to_string(index=False))

# ---- plot: paired lines for epithelial and stromal cells ----
want = [t for t in g[type_c].unique() if any(k in str(t).lower() for k in ("epithelial", "stromal"))]
if want:
    fig, axes = plt.subplots(1, len(want), figsize=(4.5 * len(want), 4.5), squeeze=False)
    for ax, t in zip(axes[0], want):
        piv = g[g[type_c] == t].pivot(index=pat_c, columns=cls_c, values="log2cpm")
        if not {"Tumor", "Normal"} <= set(piv.columns):
            ax.set_title(str(t) + " (no data)")
            continue
        for _, r in piv.iterrows():
            if pd.notna(r["Normal"]) and pd.notna(r["Tumor"]):
                ax.plot([0, 1], [r["Normal"], r["Tumor"]], color="grey", marker="o", lw=1)
            elif pd.notna(r["Tumor"]):
                ax.plot([1], [r["Tumor"]], color="firebrick", marker="o")
            elif pd.notna(r["Normal"]):
                ax.plot([0], [r["Normal"]], color="steelblue", marker="o")
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Normal", "Tumor"])
        ax.set_ylabel("%s pseudobulk log2 CPM" % GENE)
        ax.set_title(str(t))
    fig.tight_layout()
    fig.savefig("%s_sc_pseudobulk.png" % GENE, dpi=200)
print("Done. Saved %s_sc_pseudobulk_*.csv and %s_sc_pseudobulk.png" % (GENE, GENE))
