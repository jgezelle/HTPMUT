#!/usr/bin/env python3
"""
Generate sequence logos from mutation frequencies with SD shading.
Requires the output from mut_freq_norm.py.
"""

import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

# ---------------------------
# 1. Set directories and files
# ---------------------------
INPUT_FILE = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo/mutation_summary_per_pos.tsv"
OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo/logos"
os.makedirs(OUTPUT_DIR, exist_ok=True)

summary = pd.read_csv(INPUT_FILE, sep="\t")

# ---------------------------
# 2. Generate sequence logos
# ---------------------------
for virus in summary["virus"].unique():
    for condition in summary["condition"].unique():
        df_plot = summary[(summary["virus"]==virus) & (summary["condition"]==condition)].copy()
        if df_plot.empty:
            continue

        # Build matrix for logomaker (A,C,G,T frequencies)
        df_logo = pd.DataFrame(columns=["A","C","G","T"])
        for idx, row in df_plot.iterrows():
            freqs = {"A":0,"C":0,"G":0,"T":0}
            ref_base = "A"  # default in case ref column missing
            if "ref" in row:
                ref_base = row["ref"]
            freqs[ref_base] = 1 - row["mean_mut_rate"]
            mut_per_base = row["mean_mut_rate"]/3
            for b in freqs:
                if b != ref_base:
                    freqs[b] = mut_per_base
            df_logo = pd.concat([df_logo, pd.DataFrame([freqs])], ignore_index=True)

        # Make logo
        logo = logomaker.Logo(df_logo, shade_below=0.5, fade_below=0.5, font_name='Arial')
        plt.title(f"{virus} {condition}")
        plt.xlabel("Position")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(os.path.join(OUTPUT_DIR, f"{virus}_{condition}_logo.png"))
        plt.close()

print(f"Done! Logos saved in {OUTPUT_DIR}")
