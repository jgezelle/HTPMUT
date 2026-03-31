#!/usr/bin/env python3
"""
plot_mutindex.py  —  visualize per-position mutagenesis index

Reads the mutindex.tsv output from wtnorm_fc.py and produces:
  1. Bar plot  — x=position, y=mean_log2_mi, colored by depletion/enrichment
  2. Single-row heatmap  — same color scale as main heatmap, one row
  3. VARNA color file  — per-position score for superimposing on RNA structure

Mutagenesis index interpretation:
  log2_mi < 0  →  mutations depleted in Xrn1+ (position is structurally intolerant)
  log2_mi = 0  →  neutral (mutations equally represented)
  log2_mi > 0  →  mutations enriched (rare, possible artifact or compensatory)

VARNA color file format:
  One line per nucleotide: position<tab>value
  Load in VARNA with: File > Load > Annotation > Custom Colors

Usage:
    python plot_mutindex.py \\
        --tsv    scnv_mutindex.tsv \\
        --title  "scnv XR2" \\
        --out    scnv_mutindex

    # outputs:
    #   scnv_mutindex_barplot.pdf
    #   scnv_mutindex_heatmap.pdf
    #   scnv_mutindex.varna     (VARNA color annotation file)

Options:
    --vmax        color/bar scale max (symmetric, default 6.0)
    --region      START END  positions to plot (1-based)
    --p-thresh    significance threshold for marking bars (default 0.05)
    --fc-thresh   |log2_mi| threshold for marking (default 1.0)
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

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'

BASE_TILE_COLORS = {
    'A': '#008B8B',   # teal
    'C': '#800000',   # maroon
    'G': '#2E8B57',   # emerald
    'T': '#B8860B',   # gold
}


# ── load data ─────────────────────────────────────────────────────────────────

def load_tsv(path: str, region=None) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    if region is not None:
        df = df[(df['pos'] >= region[0]) & (df['pos'] <= region[1])]
    return df.sort_values('pos').reset_index(drop=True)


# ── bar plot ──────────────────────────────────────────────────────────────────

def plot_barplot(df: pd.DataFrame, title: str, out_path: str,
                 vmax: float = 6.0,
                 p_thresh: float = 0.05,
                 fc_thresh: float = 1.0):
    """
    Bar plot: x = position, y = mean_log2_mi
    Bars colored blue (depleted) or red (enriched).
    Significant bars (p_adj_mi < p_thresh AND |log2_mi| > fc_thresh)
    marked with a dot at the top.
    Error bars show replicate spread (min/max across reps).
    """
    fig, ax = plt.subplots(figsize=(max(8, len(df) * 0.22), 4), dpi=150)

    log2mi_cols = [c for c in df.columns if c.startswith('log2_mi_rep')]

    positions = df['pos'].values
    means     = df['mean_log2_mi'].values
    low_flag  = df['low_depth_flag'].values

    # error bar: min and max across replicates
    if len(log2mi_cols) > 1:
        rep_mat  = df[log2mi_cols].values
        err_low  = means - rep_mat.min(axis=1)
        err_high = rep_mat.max(axis=1) - means
        yerr     = [err_low, err_high]
    else:
        yerr = None

    # colors by direction
    colors = ['#2166AC' if v < 0 else '#B2182B' for v in means]
    # low depth → grey
    colors = ['#cccccc' if low_flag[i] else colors[i]
              for i in range(len(colors))]

    x = np.arange(len(positions))

    ax.bar(x, means,
           color=colors, width=0.8,
           yerr=yerr,
           error_kw=dict(elinewidth=0.6, ecolor='#555555', capsize=1.5),
           zorder=2)

    # significance dots
    sig = (
        (~df['low_depth_flag']) &
        (df['p_adj_mi'] < p_thresh) &
        (df['mean_log2_mi'].abs() > fc_thresh)
    )
    for i, (is_sig, val) in enumerate(zip(sig, means)):
        if is_sig:
            y_dot = val + (0.15 if val >= 0 else -0.25)
            ax.plot(x[i], y_dot, 'k.', markersize=3, zorder=3)

    # WT base row beneath x-axis
    ax.set_xticks(x)
    ax.set_xticklabels([])   # hide default labels; draw custom below

    for i, (pos, wt_base) in enumerate(zip(positions, df['wt_base'])):
        color = BASE_TILE_COLORS.get(wt_base, '#cccccc')
        ax.annotate('', xy=(x[i], -vmax * 0.08),
                    xycoords='data',
                    xytext=(x[i], -vmax * 0.18),
                    arrowprops=dict(arrowstyle='-', lw=0))
        ax.add_patch(mpatches.FancyBboxPatch(
            (x[i] - 0.4, -vmax * 0.22), 0.8, vmax * 0.14,
            boxstyle='square,pad=0',
            facecolor=color, edgecolor='white', lw=0.3,
            clip_on=False, zorder=4
        ))
        ax.text(x[i], -vmax * 0.15, wt_base,
                ha='center', va='center',
                fontsize=max(4, min(8, int(180 / len(positions)))),
                fontweight='bold', color='white',
                clip_on=False, zorder=5)
        if pos % 5 == 0 or i == 0:
            ax.text(x[i], -vmax * 0.30, str(pos),
                    ha='center', va='top',
                    fontsize=max(4, min(7, int(150 / len(positions)))),
                    color='#333333', clip_on=False, zorder=5)

    ax.axhline(0, color='#333333', lw=0.8, zorder=1)
    ax.axhline( fc_thresh, color='#aaaaaa', lw=0.6, ls='--', zorder=1)
    ax.axhline(-fc_thresh, color='#aaaaaa', lw=0.6, ls='--', zorder=1)

    ax.set_ylim(-vmax, vmax)
    ax.set_xlim(-0.5, len(positions) - 0.5)
    ax.set_ylabel('mean log₂(mutagenesis index)\n← depleted  |  enriched →',
                  fontsize=9)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=8)

    # legend
    dep_patch = mpatches.Patch(facecolor='#2166AC', label='depleted (structurally intolerant)')
    enr_patch = mpatches.Patch(facecolor='#B2182B', label='enriched')
    lo_patch  = mpatches.Patch(facecolor='#cccccc', label='low depth')
    dot_patch = mpatches.Patch(facecolor='black',
                               label=f'p_adj < {p_thresh}, |log₂MI| > {fc_thresh}')
    ax.legend(handles=[dep_patch, enr_patch, lo_patch, dot_patch],
              fontsize=7, frameon=False, loc='lower right')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight',
                dpi=200 if out_path.endswith('.png') else None)
    print(f"Saved bar plot: {out_path}")
    plt.close()


# ── single-row heatmap ────────────────────────────────────────────────────────

def plot_heatmap_row(df: pd.DataFrame, title: str, out_path: str,
                     vmax: float = 6.0):
    """
    Single-row heatmap of mutagenesis index — matches style of per-base heatmap.
    Color scale: blue = depleted, red = enriched, white = neutral.
    """
    positions = df['pos'].values
    values    = df['mean_log2_mi'].values
    wt_bases  = df['wt_base'].values
    low_flag  = df['low_depth_flag'].values
    n_pos     = len(positions)

    fig_w = max(10, n_pos * 0.25)
    fig   = plt.figure(figsize=(fig_w, 2.8), dpi=150)
    gs    = fig.add_gridspec(
        3, 1,
        height_ratios=[1.8, 0.45, 0.22],
        hspace=0.06,
        left=0.07, right=0.94,
        top=0.82, bottom=0.10
    )
    ax_heat = fig.add_subplot(gs[0])
    ax_wt   = fig.add_subplot(gs[1])
    ax_cbar = fig.add_subplot(gs[2])

    cmap = matplotlib.colormaps['RdBu_r']
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # build 1 x n_pos matrix
    mat = values.reshape(1, -1).copy()
    # grey out low depth
    rgba = np.zeros((1, n_pos, 4))
    for pi in range(n_pos):
        if low_flag[pi] or np.isnan(mat[0, pi]):
            rgba[0, pi] = [0.85, 0.85, 0.85, 0.5]
        else:
            rgba[0, pi, :3] = cmap(norm(mat[0, pi]))[:3]
            rgba[0, pi,  3] = 1.0

    ax_heat.imshow(rgba, aspect='auto', interpolation='nearest',
                   extent=[-0.5, n_pos - 0.5, 0.5, -0.5])
    ax_heat.set_yticks([0])
    ax_heat.set_yticklabels(['MI'], fontsize=9, fontweight='bold')
    ax_heat.set_xticks([])
    ax_heat.set_xlim(-0.5, n_pos - 0.5)
    ax_heat.set_title(f'{title} — mutagenesis index', fontsize=11,
                      fontweight='bold', pad=6)

    # WT base row
    ax_wt.set_xlim(-0.5, n_pos - 0.5)
    ax_wt.set_ylim(0, 1)
    ax_wt.axis('off')
    base_fs = max(5, min(11, int(200 / n_pos)))
    pos_fs  = max(4, min(8,  int(160 / n_pos)))

    for pi, (pos, wt_base) in enumerate(zip(positions, wt_bases)):
        color = BASE_TILE_COLORS.get(wt_base, '#cccccc')
        ax_wt.add_patch(mpatches.Rectangle(
            (pi - 0.5, 0.05), 1, 0.90,
            facecolor=color, edgecolor='white', lw=0.4, zorder=1))
        ax_wt.text(pi, 0.68, wt_base,
                   ha='center', va='center',
                   fontsize=base_fs, fontweight='bold',
                   color='white', zorder=2)
        if pos % 5 == 0 or pi == 0:
            ax_wt.text(pi, 0.20, str(pos),
                       ha='center', va='center',
                       fontsize=pos_fs, color='white', zorder=2)
    ax_wt.text(-0.8, 0.5, 'WT', ha='right', va='center',
               fontsize=9, fontweight='bold')

    # colorbar
    sm = matplotlib.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cb = fig.colorbar(sm, cax=ax_cbar, orientation='horizontal')
    cb.set_label(
        f'mean log₂(mutagenesis index)   |   scale: ±{vmax}   |   '
        f'blue = structurally intolerant,  red = enriched',
        fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    plt.savefig(out_path, bbox_inches='tight',
                dpi=200 if out_path.endswith('.png') else None)
    print(f"Saved heatmap row: {out_path}")
    plt.close()


# ── VARNA color file ──────────────────────────────────────────────────────────

def write_varna_colors(df: pd.DataFrame, out_path: str, vmax: float = 6.0):
    """
    Write a VARNA-compatible color annotation file.

    Format: one line per nucleotide
        position<tab>value
    where value is mean_log2_mi clamped to [-vmax, vmax].

    In VARNA: Annotations > Color map > Load custom values
    The color scale in VARNA should be set to match vmax.

    Positions with low_depth_flag are written as 0 (neutral gray in VARNA).
    """
    with open(out_path, 'w') as fh:
        fh.write(f"# VARNA color annotation — mutagenesis index\n")
        fh.write(f"# mean log2(WT-normalized mutagenesis index)\n")
        fh.write(f"# scale: -{vmax} to +{vmax}\n")
        fh.write(f"# negative = structurally intolerant (depleted mutations)\n")
        fh.write(f"# generated by plot_mutindex.py\n")
        fh.write("pos\tlog2_mut_index\n")
        for _, row in df.iterrows():
            val = 0.0 if row['low_depth_flag'] else float(
                np.clip(row['mean_log2_mi'], -vmax, vmax))
            fh.write(f"{int(row['pos'])}\t{val:.4f}\n")
    print(f"Saved VARNA file: {out_path}")


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tsv',      required=True,
                   help='mutindex.tsv from wtnorm_fc.py')
    p.add_argument('--title',    default='Mutagenesis index')
    p.add_argument('--out',      default='mutindex',
                   help='Output path stem (no extension)')
    p.add_argument('--vmax',     type=float, default=6.0)
    p.add_argument('--region',   nargs=2, type=int, metavar=('START', 'END'),
                   default=None)
    p.add_argument('--p-thresh', type=float, default=0.05)
    p.add_argument('--fc-thresh',type=float, default=1.0,
                   help='|log2_mi| threshold for significance dots (default 1.0)')
    return p.parse_args()


def main():
    args   = parse_args()
    region = tuple(args.region) if args.region else None
    ext    = '.pdf'

    print(f"Loading {args.tsv}...")
    df = load_tsv(args.tsv, region=region)
    print(f"  {len(df)} positions")

    plot_barplot(df, title=args.title,
                 out_path=f"{args.out}_barplot{ext}",
                 vmax=args.vmax,
                 p_thresh=args.p_thresh,
                 fc_thresh=args.fc_thresh)

    plot_heatmap_row(df, title=args.title,
                     out_path=f"{args.out}_heatmap{ext}",
                     vmax=args.vmax)

    write_varna_colors(df, out_path=f"{args.out}.varna",
                       vmax=args.vmax)


if __name__ == '__main__':
    main()
