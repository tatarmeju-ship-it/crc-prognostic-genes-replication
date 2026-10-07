"""
Cut the DCBLD2 region (chromosome 3, 96-102 Mb; a wide window that is valid
for both GRCh37 and GRCh38) out of a big GWAS summary-statistics file.

Usage (Anaconda Prompt):
    conda activate base
    cd /d D:\\crc
    python extract_gwas_region.py GCST90255675_buildGRCh37.tsv.gz

Use the real file name you downloaded. Output: DCBLD2_region_gwas.tsv
"""
import sys
import pandas as pd

if len(sys.argv) < 2:
    sys.exit("Give the GWAS file name, for example: python extract_gwas_region.py FILE.tsv.gz")
path = sys.argv[1]
LO, HI = 96_000_000, 102_000_000


def find(cols, *names):
    low = {c.lower().lstrip("#"): c for c in cols}
    for n in names:
        if n.lower() in low:
            return low[n.lower()]
    return None


first = pd.read_csv(path, sep="\t", nrows=3, dtype=str)
cols = list(first.columns)
print("Columns:", cols)
print(first.to_string())

chrom = find(cols, "chromosome", "chr", "chrom", "CHR")
pos = find(cols, "base_pair_location", "pos", "position", "bp", "BP")
if chrom is None or pos is None:
    sys.exit("Could not find chromosome / position columns. Send me the column list above.")

kept, total = [], 0
for chunk in pd.read_csv(path, sep="\t", chunksize=500_000, dtype=str, low_memory=False):
    total += len(chunk)
    c = chunk[chrom].astype(str).str.replace("chr", "", regex=False)
    p = pd.to_numeric(chunk[pos], errors="coerce")
    sel = chunk[(c == "3") & (p >= LO) & (p <= HI)]
    if len(sel):
        kept.append(sel)
    print("read %d rows, kept so far %d" % (total, sum(len(k) for k in kept)), flush=True)

if not kept:
    sys.exit("No variants found in the window. Send me the column list and a few rows.")
out = pd.concat(kept)
out.to_csv("DCBLD2_region_gwas.tsv", sep="\t", index=False)
print("Saved DCBLD2_region_gwas.tsv with %d variants (of %d total)." % (len(out), total))
