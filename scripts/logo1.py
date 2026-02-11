#!/usr/bin/env python3
"""
Create sequence logos showing mean mutation frequency per position,
with shading indicating SD across replicates.
Uses the output from mut_freq_norm.py
"""

import os
import pandas as pd
import logomaker
import matplotlib.pyplot as plt
import numpy as np

OUTPUT_DIR = "/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_logo"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Load summary table
summary = pd.read_csv(os.path.join(OUTPUT_DIR, "mutation_summary_per_pos.tsv"), sep="\t")

# Loop through virus & condition groups
for (virus, condition), df_group in summary.groupby(["virus","condition"]):
    
    positions = df_group["pos"].values
    # Create empty matrix for A, C, G, T
    logo_matrix = pd.DataFrame(0, index=positions, columns=["A","C","G","T"])
    
    # Fill in mean mutation frequencies for non-reference bases
    for i, row in df_group.iterrows():
        mut = row["mean_mut_rate"]
        sd = row["sd_mut_rate"]
        ref = row.get("ref", None)  # add ref to summary if needed
        
        # Assign mut freq to non-ref bases
        other_bases = [b for b in ["A","C","G","T"] if b != ref]
        for b in other_bases:
            logo_matrix.at[row["pos"], b] = mut
    
    # Convert to numpy for shading
    mean_vals = logo_matrix.to_numpy()
    sd_vals = np.tile(df_group["sd_mut_rate"].to_numpy().reshape(-1,1), (1,4))
    
    # Plot
    plt.figure(figsize=(max(10, len(positions)/2), 4))
    
    # Logo
    logo = logomaker.Logo(logo_matrix,
                          shade_below=0.5,
                          fade_below=0.5)
    
    # Add SD as translucent overlay
    for j, base in enumerate(["A","C","G","T"]):
        plt.fill_between(positions, 
                         mean_vals[:, j]-sd_vals[:, j], 
                         mean_vals[:, j]+sd_vals[:, j],
                         color='grey', alpha=0.2)
    
    plt.title(f"{virus} {condition} mutation logo")
    plt.xlabel("Position")
    plt.ylabel("Mutation frequency")
    plt.ylim(0, 1)
    plt.tight_layout()
    
    plt.savefig(os.path.join(OUTPUT_DIR, f"{virus}_{condition}_mutation_logo.png"), dpi=300)
    plt.close()

print(f"logos saved in {OUTPUT_DIR}")
