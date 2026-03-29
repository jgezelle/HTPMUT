#!/usr/bin/env python3

import glob
import os
import pandas as pd

READCOUNT_DIR = os.environ["READCOUNT_DIR"]
OUTDIR = os.environ["OUTDIR"]

BASES = ["A","C","G","T"]

def parse_sample(sample):
    # e.g. 260325-r0415-S-X-2-A3477_R1.NOLEAD
    parts = sample.split("-")
    virus_code = parts[2]      # S or B
    condition  = parts[3]      # X or R
    replicate  = parts[4]      # 1 or 2

    virus = "SCNV" if virus_code == "S" else "BENY"
    return virus, condition, replicate


def parse_readcount(path):
    sample = os.path.basename(path).replace(".readcount.txt","")
    virus, condition, rep = parse_sample(sample)

    rows = []

    with open(path) as f:
        for line in f:
            cols = line.strip().split("\t")

            contig = cols[0]
            pos = int(cols[1])
            ref = cols[2]
            depth = int(cols[3])

            counts = {b:0 for b in BASES}

            for tok in cols[4:]:
                if tok[0] in BASES and tok[1]==":":
                    counts[tok[0]] = int(tok.split(":")[1])

            for alt in BASES:
                if alt == ref:
                    continue

                feature = f"{pos}_{ref}>{alt}"

                rows.append({
                    "sample": sample,
                    "virus": virus,
                    "condition": condition,
                    "replicate": rep,
                    "feature": feature,
                    "count": counts[alt]
                })

    return pd.DataFrame(rows)


# ----------------------

files = sorted(glob.glob(f"{READCOUNT_DIR}/*.readcount.txt"))

dfs = [parse_readcount(f) for f in files]
df = pd.concat(dfs)

# counts matrix
mat = df.pivot(index="feature", columns="sample", values="count").fillna(0)
mat = mat.astype(int)

mat.to_csv(f"{OUTDIR}/alt_counts.tsv", sep="\t")

# sample metadata
meta = df[["sample","virus","condition","replicate"]].drop_duplicates()
meta = meta.set_index("sample")

meta.to_csv(f"{OUTDIR}/sample_table.tsv", sep="\t")

print("Done.")