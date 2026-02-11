#!/usr/bin/env python3
"""
Compute mutation frequency per position for each sample (per-replicate).
Output per-sample mut_rate table for downstream visualization.
"""

import os
import glob
import pandas as pd

# ---------------------------
# Directories
# ---------------------------
PILEUP_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/05_pileup"
OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------
# Load files
# ---------------------------
files = glob.glob(f"{PILEUP_DIR}/*.pileup.counts.tsv")
all_data = []

for f in files:
    basename = os.path.basename(f).replace(".pileup.counts.tsv", "")

    virus = "ZIKV" if "-Z" in basename else "DENV" if "-D" in basename else "Unknown"
    condition = "R" if "ZR" in basename or "DR" in basename else "X" if "ZX" in basename or "DX" in basename else "Unknown"
    replicate = basename.split('-')[2]  # e.g., 01, 02, 03

    df = pd.read_csv(f, sep="\t")
    df["depth"] = df[["A","C","G","T"]].sum(axis=1)
    df["ref_count"] = df.apply(lambda row: row[row["ref"]], axis=1)
    df["mut_rate"] = 1 - (df["ref_count"] / df["depth"])

    df_out = df[["pos", "mut_rate"]].copy()
    df_out["sample"] = basename
    df_out["virus"] = virus
    df_out["condition"] = condition
    df_out["replicate"] = replicate
    all_data.append(df_out)

combined = pd.concat(all_data, ignore_index=True)
combined.to_csv(os.path.join(OUTPUT_DIR, "mutation_per_replicate.tsv"), sep="\t", index=False)
print(f"Done! Output in {OUTPUT_DIR}")
