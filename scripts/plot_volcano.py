#!/usr/bin/env python3
"""
plot_volcano.py  —  standalone plotting script for wtnorm_fc.py output

Produces two figures from a wtnorm_fc TSV:
  1. Volcano plot  — mean log2FC (x) vs -log10(p_adj) (y)
       - y-axis capped at 50; saturated points shown as triangles at top
       - point size scaled by |log2FC| (effect size)
       - non-overlapping labels via adjustText
  2. MA plot       — mean log2(depth) (x) vs mean log2FC (y)
       - shows whether FC correlates with coverage (depth bias check)
       - same coloring and size scaling

Usage:
    python plot_volcano.py \\
        --tsv    scnv_wtnorm_fc.tsv \\
        --title  "scnv XR2" \\
        --out    scnv_plots.pdf

    # outputs: scnv_plots_volcano.pdf  and  scnv_plots_ma.pdf
    # (or .png if --out ends with .png)

Requirements:
    pip install adjustText
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from adjustText import adjust_text

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['axes.spines.top'] = False
matplotlib.rcParams['axes.spines.right'] = False

BASE_COLORS = {
    'A': '#E64B35',
    'C': '#4DBBD5',
    'G': '#F39B7F',
    'T': '#3C5488',
}

Y_CAP       = 50      # -log10(p) above this → shown as triangle at cap
P_THRESH    = 0.05
FC_THRESH   = 0.585   # log2(1.5)
SIZE_MIN    = 10
SIZE_MAX    = 120


# ── helpers ───────────────────────────────────────────────────────────────────

def load_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep='\t')
    df = df[df['mut_base'] != df['wt_base']].copy()   # exclude WT rows
    df['neg_log10_p'] = -np.log10(df['p_adj'].clip(lower=1e-300))
    df['saturated']   = df['neg_log10_p'] > Y_CAP
    df['y_plot']      = df['neg_log10_p'].clip(upper=Y_CAP)
    df['size']        = (
        SIZE_MIN + (SIZE_MAX - SIZE_MIN) *
        (df['mean_log2fc'].abs() / df['mean_log2fc'].abs().max())
    )
    # mean depth across replicates and conditions
    depth_cols = [c for c in df.columns if c.startswith('depth_')]
    df['mean_depth'] = df[depth_cols].mean(axis=1)
    df['log2_mean_depth'] = np.log2(df['mean_depth'].clip(lower=1))
    df['label'] = df['wt_base'] + df['pos'].astype(str) + df['mut_base']
    return df


def point_style(row):
    """marker and alpha by significance and saturation"""
    if row['saturated']:
        return '^', 0.85   # triangle up for capped points
    return 'o', 0.70


def sig_mask(df, p_thresh=P_THRESH, fc_thresh=FC_THRESH):
    return (
        (df['p_adj'] < p_thresh) &
        (df['mean_log2fc'].abs() > fc_thresh) &
        (~df['low_depth_flag'])
    )


# ── volcano ───────────────────────────────────────────────────────────────────

def plot_volcano(df: pd.DataFrame, title: str, out_path: str):
    fig, ax = plt.subplots(figsize=(7, 6), dpi=150)

    # threshold lines
    ax.axvline( FC_THRESH, color='#aaaaaa', lw=0.8, ls='--', zorder=1)
    ax.axvline(-FC_THRESH, color='#aaaaaa', lw=0.8, ls='--', zorder=1)
    ax.axhline(-np.log10(P_THRESH), color='#aaaaaa', lw=0.8, ls='--', zorder=1)
    ax.axvline(0, color='#dddddd', lw=0.6, zorder=1)
    # cap line
    ax.axhline(Y_CAP, color='#cccccc', lw=0.5, ls=':', zorder=1)
    ax.text(ax.get_xlim()[0] if ax.get_xlim()[0] < -1 else -0.5,
            Y_CAP + 0.5, f'>{Y_CAP}', fontsize=6, color='#888888', va='bottom')

    texts = []
    sig = sig_mask(df)

    for mut_base, grp in df.groupby('mut_base'):
        color = BASE_COLORS.get(mut_base, '#999999')

        # non-saturated solid points
        ns = grp[~grp['saturated'] & ~grp['low_depth_flag']]
        ax.scatter(ns['mean_log2fc'], ns['y_plot'],
                   c=color, s=ns['size'], alpha=0.70,
                   marker='o', lw=0, zorder=2)

        # saturated points as triangles at cap
        sat = grp[grp['saturated'] & ~grp['low_depth_flag']]
        ax.scatter(sat['mean_log2fc'], sat['y_plot'],
                   c=color, s=sat['size'], alpha=0.85,
                   marker='^', lw=0, zorder=3)

        # low depth as open circles
        lo = grp[grp['low_depth_flag']]
        if len(lo):
            ax.scatter(lo['mean_log2fc'], lo['y_plot'],
                       facecolors='none', edgecolors=color,
                       s=lo['size'] * 0.6, alpha=0.45,
                       marker='o', lw=0.7, zorder=2)

        # collect labels for significant hits
        grp_sig = grp[sig[grp.index] & ~grp['low_depth_flag']]
        for _, row in grp_sig.iterrows():
            texts.append(ax.text(
                row['mean_log2fc'], row['y_plot'], row['label'],
                fontsize=5.5, color='#1a1a1a',
                ha='center', va='bottom'
            ))

    # adjust text to avoid overlaps, with leader lines
    if texts:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.5),
            expand=(1.3, 1.5),
            force_text=(0.3, 0.5),
        )

    ax.set_xlabel('mean log₂(WT-normalized FC)', fontsize=11)
    ax.set_ylabel(f'−log₁₀(p adj, BH)  [capped at {Y_CAP}]', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=8)

    # legend: base colors + size scale + markers
    color_handles = [
        mpatches.Patch(facecolor=BASE_COLORS[b], label=f'→{b}')
        for b in ['A', 'C', 'G', 'T']
    ]
    size_handles = [
        mlines.Line2D([0],[0], marker='o', color='w',
                      markerfacecolor='#666', markersize=np.sqrt(s),
                      label=lbl)
        for s, lbl in [(SIZE_MIN, 'small FC'), (SIZE_MAX//2, 'med FC'),
                       (SIZE_MAX, 'large FC')]
    ]
    sat_handle = mlines.Line2D([0],[0], marker='^', color='w',
                               markerfacecolor='#666', markersize=7,
                               label=f'p > {Y_CAP} (saturated)')
    lo_handle  = mlines.Line2D([0],[0], marker='o', color='w',
                               markerfacecolor='none',
                               markeredgecolor='#666', markersize=6,
                               label='low depth')

    leg1 = ax.legend(handles=color_handles, title='mut base',
                     fontsize=7, title_fontsize=7,
                     frameon=False, loc='upper left',
                     bbox_to_anchor=(0, 1))
    ax.add_artist(leg1)
    ax.legend(handles=size_handles + [sat_handle, lo_handle],
              fontsize=7, frameon=False,
              loc='upper left', bbox_to_anchor=(0.12, 1))

    # counts
    enr = df[sig & (df['mean_log2fc'] >  FC_THRESH) & ~df['low_depth_flag']]
    dep = df[sig & (df['mean_log2fc'] < -FC_THRESH) & ~df['low_depth_flag']]
    ax.text(0.98, 0.98, f'enriched: {len(enr)}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=8, color='#c0392b')
    ax.text(0.02, 0.98, f'depleted: {len(dep)}',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=8, color='#2980b9')

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    print(f"Saved volcano: {out_path}")
    plt.close()


# ── MA plot ───────────────────────────────────────────────────────────────────

def plot_ma(df: pd.DataFrame, title: str, out_path: str):
    """
    MA plot: x = log2(mean depth), y = mean log2FC
    Useful for spotting depth-dependent FC bias.
    A well-behaved experiment shows no trend — FC should be independent of depth.
    """
    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)

    ax.axhline(0, color='#dddddd', lw=0.8, zorder=1)
    ax.axhline( FC_THRESH, color='#aaaaaa', lw=0.7, ls='--', zorder=1)
    ax.axhline(-FC_THRESH, color='#aaaaaa', lw=0.7, ls='--', zorder=1)

    sig = sig_mask(df)
    texts = []

    for mut_base, grp in df.groupby('mut_base'):
        color = BASE_COLORS.get(mut_base, '#999999')

        non_sig = grp[~sig[grp.index] & ~grp['low_depth_flag']]
        ax.scatter(non_sig['log2_mean_depth'], non_sig['mean_log2fc'],
                   c=color, s=non_sig['size'] * 0.5,
                   alpha=0.35, lw=0, marker='o', zorder=2)

        sig_grp = grp[sig[grp.index] & ~grp['low_depth_flag']]
        ax.scatter(sig_grp['log2_mean_depth'], sig_grp['mean_log2fc'],
                   c=color, s=sig_grp['size'],
                   alpha=0.85, lw=0, marker='o', zorder=3)

        lo = grp[grp['low_depth_flag']]
        if len(lo):
            ax.scatter(lo['log2_mean_depth'], lo['mean_log2fc'],
                       facecolors='none', edgecolors=color,
                       s=lo['size'] * 0.5, alpha=0.4,
                       marker='o', lw=0.7, zorder=2)

        for _, row in sig_grp.iterrows():
            texts.append(ax.text(
                row['log2_mean_depth'], row['mean_log2fc'], row['label'],
                fontsize=5.5, color='#1a1a1a', ha='center', va='bottom'
            ))

    if texts:
        adjust_text(
            texts, ax=ax,
            arrowprops=dict(arrowstyle='-', color='#999999', lw=0.5),
            expand=(1.3, 1.5),
        )

    ax.set_xlabel('log₂(mean depth)', fontsize=11)
    ax.set_ylabel('mean log₂(WT-normalized FC)', fontsize=11)
    ax.set_title(f'{title} — MA plot', fontsize=12, fontweight='bold', pad=8)

    color_handles = [
        mpatches.Patch(facecolor=BASE_COLORS[b], label=f'→{b}')
        for b in ['A', 'C', 'G', 'T']
    ]
    ax.legend(handles=color_handles, title='mut base',
              fontsize=7, title_fontsize=7,
              frameon=False, loc='upper left')

    # LOESS-style trend line to spot depth bias
    try:
        from statsmodels.nonparametric.smoothers_lowess import lowess
        valid = df[~df['low_depth_flag']].dropna(
            subset=['log2_mean_depth', 'mean_log2fc'])
        smoothed = lowess(valid['mean_log2fc'], valid['log2_mean_depth'],
                          frac=0.4, return_sorted=True)
        ax.plot(smoothed[:, 0], smoothed[:, 1],
                color='#333333', lw=1.2, ls='-', alpha=0.6,
                label='LOWESS trend', zorder=4)
        ax.legend(fontsize=7, frameon=False)
    except ImportError:
        pass

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight')
    print(f"Saved MA plot: {out_path}")
    plt.close()


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--tsv',   required=True,
                   help='wtnorm_fc TSV output file')
    p.add_argument('--title', default='Mutation FC',
                   help='Plot title')
    p.add_argument('--out',   default='plots.pdf',
                   help='Output path base (suffix _volcano / _ma added)')
    p.add_argument('--p-thresh',  type=float, default=P_THRESH)
    p.add_argument('--fc-thresh', type=float, default=1.5,
                   help='FC threshold (linear, default 1.5 -> log2=0.585)')
    p.add_argument('--y-cap', type=float, default=Y_CAP,
                   help=f'y-axis cap for volcano (default {Y_CAP})')
    return p.parse_args()


def main():
    args = parse_args()

    global Y_CAP, P_THRESH, FC_THRESH
    Y_CAP      = args.y_cap
    P_THRESH   = args.p_thresh
    FC_THRESH  = np.log2(args.fc_thresh)

    ext  = Path(args.out).suffix or '.pdf'
    stem = str(Path(args.out).with_suffix(''))

    print(f"Loading {args.tsv}...")
    df = load_tsv(args.tsv)
    print(f"  {len(df)} mutations loaded")

    plot_volcano(df, title=args.title,
                 out_path=f"{stem}_volcano{ext}")
    plot_ma(df, title=args.title,
            out_path=f"{stem}_ma{ext}")


if __name__ == '__main__':
    main()
