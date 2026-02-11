#!/usr/bin/env python3
import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt

OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MUT_FILE = os.path.join(OUTPUT_DIR, "mutation_per_replicate.tsv")

# IMPORTANT: paste FULL sequences here (no "...")
REF_SEQS = {
    "ZIKV": "GGATTAATATAATCTGGGAAACCAAGCTCATAGTCAGGCCGAGAACGCCATGGCACGGAAGAAGCCATGCTGCCTGTGAGCCCCTCAGAGGACACTGAGTCAAAAAACCCCAC",
    "DENV": "GGATTAATATAATCCAAGGACGTTAAAAGAAGTCAGGCCATCATAAATGCCATAGCTGGAGTAAACTATGCAGCCTGTAGCTCCACCTGAGAAGGTGTAAAAAATCCGGGAGG",
}

mut = pd.read_csv(MUT_FILE, sep="\t")

for virus, ref_seq in REF_SEQS.items():
    for condition in ["R","X"]:
        sub = mut[(mut["virus"] == virus) & (mut["condition"] == condition)].copy()
        if sub.empty:
            continue

        # mean ± SD across replicates at each position
        stats = sub.groupby("pos")["mut_rate"].agg(["mean","std","count"]).reset_index()

        # sanity check: reference length must cover max pos
        maxpos = int(stats["pos"].max())
        if len(ref_seq) < maxpos:
            raise ValueError(f"{virus} reference too short: len={len(ref_seq)} but max pos in data={maxpos}")

        # --- Build a STANDARD logo from the reference (one-hot) ---
        bases = list(ref_seq[:maxpos])  # match plotted length to data
        df_logo = logomaker.alignment_to_matrix([ref_seq[:maxpos]], to_type="counts")
        # counts from 1 sequence → one-hot per position

        plt.figure(figsize=(max(12, maxpos/4), 4))

        # Logo (WT)
        logomaker.Logo(df_logo)

        # Overlay mutation rate (mean ± SD)
        x = stats["pos"].values
        y = stats["mean"].values
        ysd = stats["std"].fillna(0).values

        plt.plot(x, y, linewidth=1.5, label="mean mutation rate")
        plt.fill_between(x, y-ysd, y+ysd, alpha=0.2, label="± SD")

        plt.title(f"{virus} {condition}: WT logo + mutation rate")
        plt.xlabel("Position")
        plt.ylabel("WT counts (logo) / mutation rate (line)")
        plt.legend()
        plt.tight_layout()

        outpng = os.path.join(OUTPUT_DIR, f"{virus}_{condition}_WTlogo_mutrate.png")
        plt.savefig(outpng, dpi=300)
        plt.close()

        print("Saved:", outpng)
