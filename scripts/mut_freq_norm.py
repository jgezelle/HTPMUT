#!/usr/bin/env python3
"""
Compute mutation frequency per position for each virus & condition.
Combines replicates and outputs mean ± SD per position.
"""

import os
import glob
import pandas as pd

# ---------------------------
# 1. Set directories
# ---------------------------
PILEUP_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/05_pileup"
OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------------------
# 2. Load files
# ---------------------------
files = glob.glob(f"{PILEUP_DIR}/*.pileup.counts.tsv")

all_data = []

for f in files:
    # Extract sample name
    basename = os.path.basename(f).replace(".pileup.counts.tsv", "")
    
    # Determine virus and condition from name
    if "Z" in basename:
        virus = "ZIKV"
    elif "D" in basename:
        virus = "DENV"
    else:
        virus = "Unknown"

    if "R" in basename:
        condition = "R"  # input / no exonuclease
    elif "X" in basename:
        condition = "X"  # exonuclease treated
    else:
        condition = "Unknown"
    
    # Load counts table
    df = pd.read_csv(f, sep="\t")
    
    # Compute depth and mutation frequency
    df["depth"] = df[["A", "C", "G", "T"]].sum(axis=1)
    
    # ---------------------------
    # Modern replacement for deprecated df.lookup
    # ---------------------------
    df["ref_count"] = df.apply(lambda row: row[row["ref"]], axis=1)
    
    # mutation rate = 1 - fraction of reads matching ref
    df["mut_rate"] = 1 - (df["ref_count"] / df["depth"])
    
    # store relevant info
    df_out = df[["pos", "mut_rate"]].copy()
    df_out["sample"] = basename
    df_out["virus"] = virus
    df_out["condition"] = condition
    
    all_data.append(df_out)

# Combine all samples into one dataframe
combined = pd.concat(all_data, ignore_index=True)

# ---------------------------
# 3. Compute mean & SD per position
# ---------------------------
summary = combined.groupby(["virus", "condition", "pos"])["mut_rate"].agg(
    mean_mut_rate="mean",
    sd_mut_rate="std",
    n_replicates="count"
).reset_index()

# ---------------------------
# 4. Write outputs
# ---------------------------
summary.to_csv(os.path.join(OUTPUT_DIR, "mutation_summary_per_pos.tsv"), sep="\t", index=False)
combined.to_csv(os.path.join(OUTPUT_DIR, "mutation_all_replicates_per_pos.tsv"), sep="\t", index=False)

print(f"Done! Outputs in {OUTPUT_DIR}")
