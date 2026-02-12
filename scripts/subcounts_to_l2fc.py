#!/usr/bin/env python3
# subcounts_to_l2fc.py
#
# Compute substitution log2 fold-changes (X vs R) per virus from substitution_counts.tsv.
# Uses frequency (alt_count / depth) per replicate, then averages within condition.
#
# Output: substitution_l2fc.tsv with columns:
#   virus, feature, pos, ref, alt, mean_freq_R, mean_freq_X, log2fc

import os
import argparse
import math
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_tsv", required=True, help="substitution_counts.tsv")
    ap.add_argument("--out_tsv", required=True, help="substitution_l2fc.tsv")
    ap.add_argument("--alpha", type=float, default=1e-6, help="Pseudocount on frequency")
    ap.add_argument("--min_mean_depth", type=float, default=0.0,
                    help="Optional: drop features with mean depth below this (after grouping)")
    args = ap.parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t")

    # keep only known viruses/conditions
    df = df[df["virus"].isin(["ZIKV", "DENV"])].copy()
    df = df[df["condition"].isin(["R", "X"])].copy()

    # frequency per row (replicate)
    df["freq"] = df["alt_count"] / df["depth"]

    # mean across replicates per condition
    g = (df.groupby(["virus", "condition", "feature", "pos", "ref", "alt"], as_index=False)
           .agg(mean_freq=("freq", "mean"),
                mean_depth=("depth", "mean"),
                mean_alt_count=("alt_count", "mean")))

    # optional depth filter
    if args.min_mean_depth > 0:
        g = g[g["mean_depth"] >= args.min_mean_depth].copy()

    # pivot to wide: columns R and X
    piv = g.pivot_table(
        index=["virus", "feature", "pos", "ref", "alt"],
        columns="condition",
        values="mean_freq",
        aggfunc="first"
    ).reset_index()

    piv["R"] = piv["R"].fillna(0.0) if "R" in piv.columns else 0.0
    piv["X"] = piv["X"].fillna(0.0) if "X" in piv.columns else 0.0

    a = args.alpha
    piv["log2fc"] = piv.apply(lambda r: math.log2((r["X"] + a) / (r["R"] + a)), axis=1)

    # rename for clarity
    piv = piv.rename(columns={"R": "mean_freq_R", "X": "mean_freq_X"})

    out_dir = os.path.dirname(args.out_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    
    piv.to_csv(args.out_tsv, sep="\t", index=False)


if __name__ == "__main__":
    main()