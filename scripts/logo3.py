#!/usr/bin/env python3
"""
Create sequence logos and overlay mutation rates for R and X conditions.
"""

import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

REF_SEQS = {
    "ZIKV": "GGA..." ,  # put full ZIKV reference here
    "DENV": "GGA..." ,  # full DENV reference
}

# Load per-replicate mutation rates
mut_df = pd.read_csv(os.path.join(OUTPUT_DIR, "mutation_per_replicate.tsv"), sep="\t")

for virus, ref_seq in REF_SEQS.items():
    for condition in ["R", "X"]:
        # subset mutation rates
        df_subset = mut_df[(mut_df["virus"]==virus) & (mut_df["condition"]==condition)]
        if df_subset.empty:
            continue

        # compute mean ± SD across replicates
        summary = df_subset.groupby("pos")["mut_rate"].agg(["mean","std"]).reset_index()

        # Create one-hot reference matrix
        df_logo = pd.DataFrame([list(base) for base in ref_seq])
        df_logo = pd.get_dummies(df_logo.stack()).groupby(level=0).sum()

        # Make logo
        logo = logomaker.Logo(df_logo)
        plt.title(f"{virus} {condition}")
        plt.xlabel("Position")
        plt.ylabel("Base")

        # Overlay mutation rate shading
        plt.fill_between(
            summary["pos"]-1,   # matplotlib is 0-indexed
            0, 
            summary["mean"], 
            color='red', alpha=0.3, label="mean mut rate"
        )
        plt.errorbar(
            summary["pos"]-1,
            summary["mean"],
            yerr=summary["std"],
            fmt='none',
            ecolor='darkred',
            alpha=0.6,
            label="SD"
        )

        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{virus}_{condition}_logo_mut.png"))
        plt.close()
