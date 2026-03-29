#!/usr/bin/env python3
"""
plot_heatmap.py  —  per-nucleotide mutation enrichment heatmap

Reads from wtnorm_fc.py TSV output.

Layout (top to bottom):
  [optional] dot-bracket structure arc diagram
  heatmap: x=position, y=ACGT, color=mean_log2(WT-normalized FC)
  WT base row: colored tile per position showing WT nucleotide identity
  colorbar

Color scale: RdBu_r diverging colormap centered at 0
  Red  = enriched in treated (Xrn1+) vs control, relative to WT
  Blue = depleted in treated vs control, relative to WT
  White = no change relative to WT (FC ~ 1)

WT cell in heatmap outlined in black (FC = 1 by construction, not shown).

Usage:
    python plot_heatmap.py \\
        --tsv    scnv_wtnorm_fc.tsv \\
        --title  "scnv XR2" \\
        --out    scnv_heatmap.pdf

    # with dot-bracket structure:
    python plot_heatmap.py \\
        --tsv       scnv_wtnorm_fc.tsv \\
        --structure "(((...)))" \\
        --title     "scnv XR2" \\
        --out       scnv_heatmap.pdf

    # zoom to region:
    python plot_heatmap.py \\
        --tsv    scnv_wtnorm_fc.tsv \\
        --region 30 80 \\
        --out    scnv_heatmap_core.pdf

    # adjust color scale:
    python plot_heatmap.py --tsv scnv_wtnorm_fc.tsv --vmax 6

Statistics shown:
    Color = mean log2(WT-normalized FC) across replicate pairs.
    FC(B,i) = [freq(B,X)/freq(B,R)] x [freq(WT,R)/freq(WT,X)]
    where X = Xrn1-treated (NOLEAD), R = control (HASLEAD).
    WT base at each position normalizes to FC=1 by construction.
    Replicates averaged after independent per-replicate FC calculation.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.ticker as ticker
from matplotlib.colors import TwoSlopeNorm
from matplotlib.path import Path as MPath

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'

BASES = ['A', 'C', 'G', 'T']

# colors for WT base row tiles
# A=Teal, C=Maroon, G=Emerald, T=Gold
BASE_TILE_COLORS = {
    'A': '#008B8B',   # teal
    'C': '#800000',   # maroon
    'G': '#2E8B57',   # emerald green
    'T': '#B8860B',   # dark gold
}


# ── data loading ──────────────────────────────────────────────────────────────

def load_tsv(path: str, region=None) -> tuple:
    """
    Load wtnorm_fc TSV.
    Returns:
        mat       np.ndarray (4, n_pos)  mean_log2fc  (NaN for WT cell)
        wt_bases  list of str            WT base at each position
        positions np.ndarray             genomic positions (1-based)
    """
    df = pd.read_csv(path, sep='\t')
    # exclude WT rows (trivially FC=1, not informative for heatmap color)
    df_mut = df[df['mut_base'] != df['wt_base']].copy()

    if region is not None:
        df_mut = df_mut[(df_mut['pos'] >= region[0]) &
                        (df_mut['pos'] <= region[1])]

    positions = np.array(sorted(df_mut['pos'].unique()))
    n_pos     = len(positions)

    mat      = np.full((4, n_pos), np.nan)
    wt_bases = []

    for pi, pos in enumerate(positions):
        sub     = df_mut[df_mut['pos'] == pos]
        wt_base = sub['wt_base'].iloc[0]
        wt_bases.append(wt_base)

        for bi, base in enumerate(BASES):
            row = sub[sub['mut_base'] == base]
            if len(row) == 0:
                continue
            mat[bi, pi] = row['mean_log2fc'].values[0]

    return mat, wt_bases, positions


# ── structure drawing ─────────────────────────────────────────────────────────

def draw_structure(ax, structure: str, n_pos: int):
    ax.set_xlim(-0.5, n_pos - 0.5)
    ax.set_ylim(0, 1)
    ax.axis('off')

    struct = structure[:n_pos]
    pairs  = {}
    stack  = []
    for i, c in enumerate(struct):
        if c == '(':
            stack.append(i)
        elif c == ')' and stack:
            j = stack.pop()
            pairs[j] = i
            pairs[i] = j

    drawn  = set()
    y_base = 0.15

    for i, c in enumerate(struct):
        if i >= n_pos:
            break
        if c == '.':
            ax.plot([i, i], [y_base, y_base - 0.10],
                    color='#555555', lw=0.9, solid_capstyle='round')
        elif c in '()' and i in pairs:
            j = pairs[i]
            if (i, j) not in drawn and j < n_pos:
                drawn.add((i, j))
                drawn.add((j, i))
                span  = abs(j - i)
                arc_h = min(0.70, 0.08 + span * 0.020)
                verts = [(i, y_base), (i, y_base + arc_h),
                         (j, y_base + arc_h), (j, y_base)]
                codes = [MPath.MOVETO, MPath.CURVE4,
                         MPath.CURVE4,  MPath.CURVE4]
                path  = MPath(verts, codes)
                ax.add_patch(mpatches.PathPatch(
                    path, fc='none', ec='#222222', lw=1.0))


# ── main plot ─────────────────────────────────────────────────────────────────

def plot_heatmap(mat, wt_bases, positions,
                 structure, title, out_path, vmax=4.0):

    n_pos      = len(positions)
    has_struct = bool(structure)

    # ── figure layout ─────────────────────────────────────────────────────────
    # rows: [structure?] | heatmap | WT base row | colorbar
    fig_w = max(10, n_pos * 0.25)
    fig_h = 4.2 + (0.7 if has_struct else 0)

    hr = []
    if has_struct:
        hr.append(0.55)   # structure
    hr += [3.2, 0.45, 0.22]   # heatmap | WT row | colorbar

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)
    gs  = fig.add_gridspec(
        len(hr), 1,
        height_ratios=hr,
        hspace=0.06,
        left=0.07, right=0.94,
        top=0.89, bottom=0.08
    )

    row_idx = 0
    if has_struct:
        ax_struct = fig.add_subplot(gs[row_idx]); row_idx += 1
    ax_heat  = fig.add_subplot(gs[row_idx]); row_idx += 1
    ax_wt    = fig.add_subplot(gs[row_idx]); row_idx += 1
    ax_cbar  = fig.add_subplot(gs[row_idx])

    # ── colormap ──────────────────────────────────────────────────────────────
    cmap = matplotlib.colormaps['RdBu_r']
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # ── heatmap image ─────────────────────────────────────────────────────────
    ax_heat.imshow(
        mat,
        aspect='auto',
        cmap=cmap,
        norm=norm,
        interpolation='nearest',
        extent=[-0.5, n_pos - 0.5, 3.5, -0.5]
    )

    # WT cell outlines (black border, no fill)
    for pi, wt_base in enumerate(wt_bases):
        if wt_base not in BASES:
            continue
        bi = BASES.index(wt_base)
        ax_heat.add_patch(mpatches.Rectangle(
            (pi - 0.5, bi - 0.5), 1, 1,
            linewidth=1.8, edgecolor='black',
            facecolor='none', zorder=3
        ))

    # thin white separators between base rows
    for y in [0.5, 1.5, 2.5]:
        ax_heat.axhline(y, color='white', lw=0.5, zorder=2)

    ax_heat.set_yticks([0, 1, 2, 3])
    ax_heat.set_yticklabels(BASES, fontsize=10,
                             fontweight='bold')
    ax_heat.set_ylabel('mutation\nto base', fontsize=9, labelpad=4)
    ax_heat.set_xlim(-0.5, n_pos - 0.5)
    ax_heat.set_xticks([])   # x labels go on WT row below
    ax_heat.set_title(title, fontsize=13, fontweight='bold', pad=8)

    # ── WT base row ───────────────────────────────────────────────────────────
    # Each column = one position; tile color = base identity; letter inside
    ax_wt.set_xlim(-0.5, n_pos - 0.5)
    ax_wt.set_ylim(0, 1)
    ax_wt.axis('off')

    for pi, (pos, wt_base) in enumerate(zip(positions, wt_bases)):
        color = BASE_TILE_COLORS.get(wt_base, '#cccccc')
        ax_wt.add_patch(mpatches.Rectangle(
            (pi - 0.5, 0.05), 1, 0.90,
            facecolor=color, edgecolor='white',
            linewidth=0.4, zorder=1
        ))
        # WT base letter — sized to fill the tile
        base_fs = max(7, min(13, int(240 / n_pos)))
        ax_wt.text(pi, 0.68, wt_base,
                   ha='center', va='center',
                   fontsize=base_fs,
                   fontweight='bold',
                   color='white', zorder=2)
        # position number — every 5th position
        if (pi == 0) or (pos % 5 == 0):
            pos_fs = max(5, min(10, int(200 / n_pos)))
            ax_wt.text(pi, 0.20, str(pos),
                       ha='center', va='center',
                       fontsize=pos_fs,
                       color='white', zorder=2)

    # row label
    ax_wt.text(-0.8, 0.5, 'WT',
               ha='right', va='center',
               fontsize=9, fontweight='bold',
               transform=ax_wt.transData)

    # ── structure ─────────────────────────────────────────────────────────────
    if has_struct:
        draw_structure(ax_struct, structure, n_pos)
        ax_struct.set_xlim(-0.5, n_pos - 0.5)

    # ── colorbar ──────────────────────────────────────────────────────────────
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cbar, orientation='horizontal')
    cb.set_label(
        f'mean log₂(WT-normalized FC)   |   scale: ±{vmax}   |   '
        f'red = enriched in Xrn1+,  blue = depleted',
        fontsize=8)
    cb.ax.tick_params(labelsize=8)
    cb.ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%.1f'))

    # ── WT outline legend ─────────────────────────────────────────────────────
    wt_patch = mpatches.Patch(facecolor='none', edgecolor='black',
                               linewidth=1.8,
                               label='WT base (FC = 1 by construction)')
    fig.legend(handles=[wt_patch],
               loc='lower right',
               bbox_to_anchor=(0.94, 0.00),
               fontsize=7.5, frameon=False)

    plt.savefig(out_path, bbox_inches='tight',
                dpi=200 if out_path.endswith('.png') else None)
    print(f"Saved: {out_path}")
    plt.close()


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tsv',       required=True,
                   help='wtnorm_fc TSV output')
    p.add_argument('--title',     default='Mutation enrichment heatmap')
    p.add_argument('--out',       default='heatmap.pdf')
    p.add_argument('--structure', default='',
                   help='Dot-bracket string (optional)')
    p.add_argument('--region',    nargs=2, type=int, metavar=('START', 'END'),
                   default=None,
                   help='Positions to plot, 1-based inclusive')
    p.add_argument('--vmax',      type=float, default=4.0,
                   help='Color scale ±max (default 4.0). Values beyond are '
                        'saturated at max color.')
    return p.parse_args()


def main():
    args   = parse_args()
    region = tuple(args.region) if args.region else None

    print(f"Loading {args.tsv}...")
    mat, wt_bases, positions = load_tsv(args.tsv, region=region)
    print(f"  {len(positions)} positions")

    plot_heatmap(
        mat=mat, wt_bases=wt_bases, positions=positions,
        structure=args.structure,
        title=args.title,
        out_path=args.out,
        vmax=args.vmax,
    )


if __name__ == '__main__':
    main()
