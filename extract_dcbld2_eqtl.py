"""
Extract the DCBLD2 rows from the eQTLGen cis-eQTL file (and the matching
allele-frequency rows) so that the small outputs can be shared and used for MR.

Put these files in the working folder (D:\\crc):
    cis-eQTL_significant_20181017.txt.gz
    2018-07-18_SNP_AF_for_AlleleB_combined_allele_c....txt.gz  (name is matched by pattern)

Run (Anaconda Prompt):   conda activate base
                         cd /d D:\\crc
                         python extract_dcbld2_eqtl.py

Outputs: DCBLD2_eqtlgen_cis.tsv  and  DCBLD2_eqtlgen_snp_af.tsv
"""
import glob
import gzip
import sys

GENE = "DCBLD2"
EQTL = "cis-eQTL_significant_20181017.txt.gz"
af_files = glob.glob("2018-07-18_SNP_AF*")


def col_index(header, *names):
    low = [h.strip().lower() for h in header]
    for n in names:
        if n.lower() in low:
            return low.index(n.lower())
    return None


# ---- 1. DCBLD2 cis-eQTLs -----------------------------------------------------
try:
    f = gzip.open(EQTL, "rt")
except OSError:
    sys.exit("Cannot open %s. Is it in this folder?" % EQTL)

header = next(f).rstrip("\n").split("\t")
print("eQTL columns:", header)
g_i = col_index(header, "GeneSymbol", "Gene Symbol", "Gene_symbol")
s_i = col_index(header, "SNP", "SNPName", "rsid")
if g_i is None or s_i is None:
    sys.exit("Could not find the gene-symbol / SNP columns. Send me the column list above.")

kept, snps = [], set()
for line in f:
    parts = line.rstrip("\n").split("\t")
    if len(parts) > g_i and parts[g_i] == GENE:
        kept.append(line)
        snps.add(parts[s_i])
f.close()
print("%s rows: %d, distinct SNPs: %d" % (GENE, len(kept), len(snps)))
if not kept:
    sys.exit("No rows for %s were found." % GENE)
with open("DCBLD2_eqtlgen_cis.tsv", "w") as out:
    out.write("\t".join(header) + "\n")
    out.writelines(kept)

# ---- 2. allele frequencies for those SNPs -------------------------------------
if not af_files:
    print("Allele-frequency file (2018-07-18_SNP_AF*) not found; skipping step 2.")
    sys.exit(0)
fa = gzip.open(af_files[0], "rt")
ah = next(fa).rstrip("\n").split("\t")
print("AF file:", af_files[0])
print("AF columns:", ah)
a_i = col_index(ah, "SNP", "SNPName", "rsid")
if a_i is None:
    a_i = 0
    print("No 'SNP' column found; assuming the first column holds SNP ids.")
n = 0
with open("DCBLD2_eqtlgen_snp_af.tsv", "w") as out:
    out.write("\t".join(ah) + "\n")
    for line in fa:
        parts = line.split("\t", a_i + 1)
        if len(parts) > a_i and parts[a_i] in snps:
            out.write(line)
            n += 1
fa.close()
print("Matched %d allele-frequency rows. Done." % n)
