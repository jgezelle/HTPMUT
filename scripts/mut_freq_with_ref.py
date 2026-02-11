#!/usr/bin/env python3
import os, glob
import pandas as pd

PILEUP_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/05_pileup"
OUT_TSV    = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo/mutation_per_replicate.tsv"
os.makedirs(os.path.dirname(OUT_TSV), exist_ok=True)

files = glob.glob(f"{PILEUP_DIR}/*.pileup.counts.tsv")
rows = []

for f in files:
    sample = os.path.basename(f).replace(".pileup.counts.tsv", "")
    parts = sample.split("-")

    # find ZR/ZX/DR/DX
    code = None
    for p in parts:
        if p in {"ZR","ZX","DR","DX"}:
            code = p
            break

    if code is None:
        virus = "Unknown"
        condition = "Unknown"
        replicate = "Unknown"
    else:
        virus = "ZIKV" if code[0] == "Z" else "DENV"
        condition = code[1]  # R or X
        idx = parts.index(code)
        replicate = parts[idx + 1] if idx + 1 < len(parts) else "Unknown"

    df = pd.read_csv(f, sep="\t")  # pos, ref, A, C, G, T
    df["depth"] = df[["A","C","G","T"]].sum(axis=1)
    df["ref_count"] = df.apply(lambda r: r[r["ref"]], axis=1)
    df["mut_rate"] = 1 - (df["ref_count"] / df["depth"])

    out = df[["pos","ref","mut_rate"]].copy()   # <-- WT base included here
    out["sample"] = sample
    out["virus"] = virus
    out["condition"] = condition
    out["replicate"] = replicate
    rows.append(out)

final = pd.concat(rows, ignore_index=True)
final.to_csv(OUT_TSV, sep="\t", index=False)

print("Wrote:", OUT_TSV)
print("Columns:", list(final.columns))
