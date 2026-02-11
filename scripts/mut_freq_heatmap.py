#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

IN_TSV = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo/mutation_per_replicate.tsv"
OUTDIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo/figs"
os.makedirs(OUTDIR, exist_ok=True)

df = pd.read_csv(IN_TSV, sep="\t")

# ensure sorting is numeric
df["pos"] = df["pos"].astype(int)
df["replicate"] = df["replicate"].astype(int)

for virus in sorted(df["virus"].unique()):
    if virus == "Unknown":
        continue

    sub = df[df["virus"] == virus].copy()

    # WT labels (pos:ref) — should be identical across samples for a virus
    wt = sub[["pos", "ref"]].drop_duplicates().sort_values("pos")
    positions = wt["pos"].values
    xlabels = [f"{p}:{b}" for p, b in zip(wt["pos"], wt["ref"])]

    # build row order: R1 R2 R3 X1 X2 X3
    row_keys = []
    for cond in ["R", "X"]:
        for rep in sorted(sub["replicate"].unique()):
            row_keys.append((cond, rep))

    # matrix: rows = (cond,rep), cols = pos
    mat = np.full((len(row_keys), len(positions)), np.nan)

    for i, (cond, rep) in enumerate(row_keys):
        r = sub[(sub["condition"] == cond) & (sub["replicate"] == rep)].sort_values("pos")
        # align by pos
        r = pd.merge(wt[["pos"]], r[["pos", "mut_rate"]], on="pos", how="left")
        mat[i, :] = r["mut_rate"].values

    # plot
    plt.figure(figsize=(max(12, len(positions)/4), 3.5))
    plt.imshow(mat, aspect="auto", interpolation="nearest")

    plt.yticks(
        ticks=np.arange(len(row_keys)),
        labels=[f"{c}{rep}" for c, rep in row_keys]
    )

    stride = 10  # change to 5 if you want denser labeling
    plt.xticks(
        ticks=np.arange(0, len(positions), stride),
        labels=[xlabels[j] for j in range(0, len(positions), stride)],
        rotation=90
    )

    plt.title(f"{virus} — mutation rate heatmap (rows=replicates, cols=pos:WT)")
    plt.xlabel("Position:WTbase")
    plt.ylabel("Condition+Replicate")
    cbar = plt.colorbar()
    cbar.set_label("Mutation rate (1 - ref fraction)")
    plt.tight_layout()

    outpng = os.path.join(OUTDIR, f"{virus}_mutrate_heatmap_reps_withWT.png")
    plt.savefig(outpng, dpi=300)
    plt.close()

    print("Saved:", outpng)
