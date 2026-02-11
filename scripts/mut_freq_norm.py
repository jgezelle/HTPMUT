#!/usr/bin/env python3
import os, glob
import pandas as pd

PILEUP_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/05_pileup"
OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

files = glob.glob(f"{PILEUP_DIR}/*.pileup.counts.tsv")
all_data = []

for f in files:
    basename = os.path.basename(f).replace(".pileup.counts.tsv", "")
    parts = basename.split("-")

    # find code among the dash-separated fields
    code = None
    for p in parts:
        if p in {"ZR","ZX","DR","DX"}:
            code = p
            break
    if code is None:
        virus = "Unknown"
        condition = "Unknown"
    else:
        virus = "ZIKV" if code[0] == "Z" else "DENV"
        condition = code[1]  # R or X

    # replicate is the token after ZR/ZX/DR/DX, e.g. ...-DR-03-...
    replicate = "Unknown"
    if code in parts:
        idx = parts.index(code)
        if idx + 1 < len(parts):
            replicate = parts[idx + 1]  # "01"/"02"/"03"

    df = pd.read_csv(f, sep="\t")
    df["depth"] = df[["A","C","G","T"]].sum(axis=1)
    df["ref_count"] = df.apply(lambda row: row[row["ref"]], axis=1)
    df["mut_rate"] = 1 - (df["ref_count"] / df["depth"])

    out = df[["pos","mut_rate"]].copy()
    out["sample"] = basename
    out["virus"] = virus
    out["condition"] = condition
    out["replicate"] = replicate
    all_data.append(out)

combined = pd.concat(all_data, ignore_index=True)
combined.to_csv(os.path.join(OUTPUT_DIR, "mutation_per_replicate.tsv"), sep="\t", index=False)

print("Wrote:", os.path.join(OUTPUT_DIR, "mutation_per_replicate.tsv"))
print("pos range:", combined["pos"].min(), "to", combined["pos"].max())
