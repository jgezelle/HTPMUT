#!/usr/bin/env python3
"""
wtnorm_fc.py  (reads directly from bam-readcount output)

WT-normalized fold change per mutation per position, across paired replicates.
Also computes a per-position mutagenesis index (aggregate mutation tolerance).

Formula (per replicate pair, per non-WT base B at position i):
    FC(B,i) = [freq(B,X) / freq(B,R)] x [freq(WT,R) / freq(WT,X)]

Mutagenesis index (per position i, per replicate):
    mut_index(i) = [1 - freq(WT,X)] / [1 - freq(WT,R)]  x  [freq(WT,R) / freq(WT,X)]
    = total non-WT frequency in X vs R, WT-normalized
    log2(mut_index) < 0 means mutations are depleted at this position (structurally intolerant)
    log2(mut_index) = 0 means neutral (mutations equally represented in X and R)

Statistics:
    Per-base FC: 2-proportion z-test on raw counts (X vs R), Stouffer-combined
                 across replicates, BH-corrected.
    Mutagenesis index: same z-test on total non-WT counts vs depth,
                       Stouffer-combined, BH-corrected separately.

Usage:
    python wtnorm_fc.py \\
        --readcount-dir  /path/to/05_readcount/260327 \\
        --virus          S \\
        --out-tsv        scnv_wtnorm_fc.tsv \\
        --out-plot       scnv_volcano.pdf \\
        --title          "scnv XR2"

Readcount filename convention:
    260325-r0415-S-X-2-A3477_R1.NOLEAD.readcount.txt
"""

import argparse
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

BASES            = ['A', 'C', 'G', 'T']
PSEUDOCOUNT      = 1
EXCLUDE_START    = 29       # skip positions 1-29 (leader + pre-structure)
MIN_DEPTH_X      = 500
MIN_DEPTH_R      = 1000
MIN_DEPTH_RATIO  = 0.5      # minimum X/R depth ratio
BASE_COLORS      = {'A': '#E64B35', 'C': '#4DBBD5', 'G': '#F39B7F', 'T': '#3C5488'}


# -- readcount parser ----------------------------------------------------------

def parse_readcount_file(path: str) -> pd.DataFrame:
    """
    Parse one bam-readcount file.
    Returns DataFrame indexed by pos with columns:
        ref  depth  A  C  G  T  freq_A  freq_C  freq_G  freq_T
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 5:
                continue
            pos   = int(cols[1])
            ref   = cols[2].upper()
            depth = int(cols[3])
            if ref not in BASES:
                continue
            counts = {b: 0 for b in BASES}
            for tok in cols[4:]:
                if not tok:
                    continue
                base = tok[0].upper()
                if base in BASES:
                    try:
                        counts[base] = int(tok.split(':')[1])
                    except (IndexError, ValueError):
                        pass
            rows.append({'pos': pos, 'ref': ref, 'depth': depth, **counts})

    df = pd.DataFrame(rows).set_index('pos')
    for b in BASES:
        df[f'freq_{b}'] = (df[b] + PSEUDOCOUNT) / (df['depth'] + 4 * PSEUDOCOUNT)
    return df


def find_sample_files(readcount_dir: str, virus: str) -> dict:
    """
    Scan readcount_dir for files matching virus code (S or B).
    Returns dict: {replicate: {'X': path, 'R': path}}
    """
    pattern = os.path.join(readcount_dir, f"*-{virus}-*.readcount.txt")
    files   = sorted(glob.glob(pattern))
    if not files:
        sys.exit(f"ERROR: No readcount files found for virus '{virus}' in {readcount_dir}")

    by_rep = defaultdict(dict)
    for f in files:
        stem  = Path(f).name.replace('.readcount.txt', '')
        parts = stem.split('-')
        cond  = parts[3]
        rep   = parts[4]
        by_rep[rep][cond] = f

    for rep, pair in list(by_rep.items()):
        if 'X' not in pair or 'R' not in pair:
            missing = 'X' if 'X' not in pair else 'R'
            print(f"WARNING: rep {rep} missing {missing} file — skipping")
            del by_rep[rep]

    if not by_rep:
        sys.exit("ERROR: No complete X/R replicate pairs found.")

    for rep, pair in sorted(by_rep.items()):
        print(f"  rep{rep}: X = {Path(pair['X']).name}")
        print(f"           R = {Path(pair['R']).name}")

    return dict(by_rep)


# -- fold change ---------------------------------------------------------------

def wt_normalized_fc(x_df: pd.DataFrame,
                     r_df: pd.DataFrame) -> tuple:
    """
    Compute per-base WT-normalized FC and per-position mutagenesis index.

    Returns:
        fc_df   — long DataFrame, one row per (pos, mut_base)
        mi_df   — wide DataFrame, one row per pos (mutagenesis index)
    """
    shared = x_df.index.intersection(r_df.index)
    fc_rows = []   # per-base results
    mi_rows = []   # per-position mutagenesis index

    for pos in shared:

        # -- position filters --------------------------------------------------
        if pos <= EXCLUDE_START:
            continue

        wt_base = x_df.loc[pos, 'ref']
        if wt_base not in BASES:
            continue

        x = x_df.loc[pos]
        r = r_df.loc[pos]

        depth_X     = int(x['depth'])
        depth_R     = int(r['depth'])
        depth_ratio = depth_X / depth_R if depth_R > 0 else 0

        if depth_ratio < MIN_DEPTH_RATIO:
            continue   # hard exclude transition zone

        low_depth = (depth_X < MIN_DEPTH_X) or (depth_R < MIN_DEPTH_R)

        freq_wt_X = x[f'freq_{wt_base}']
        freq_wt_R = r[f'freq_{wt_base}']
        wt_ratio  = freq_wt_R / freq_wt_X   # WT(R)/WT(X) normalization term

        # -- per-base FC -------------------------------------------------------
        for mut_base in BASES:
            if mut_base == wt_base:
                continue

            count_X = int(x[mut_base])
            count_R = int(r[mut_base])

            fc     = (x[f'freq_{mut_base}'] / r[f'freq_{mut_base}']) * wt_ratio
            log2fc = np.log2(fc)

            try:
                _, p_raw = proportions_ztest(
                    count=[count_X + PSEUDOCOUNT, count_R + PSEUDOCOUNT],
                    nobs =[depth_X + 4 * PSEUDOCOUNT,
                           depth_R + 4 * PSEUDOCOUNT],
                    alternative='two-sided'
                )
            except Exception:
                p_raw = np.nan

            fc_rows.append({
                'pos':         pos,
                'wt_base':     wt_base,
                'mut_base':    mut_base,
                'fc':          fc,
                'log2fc':      log2fc,
                'count_mut_X': count_X,
                'depth_X':     depth_X,
                'count_mut_R': count_R,
                'depth_R':     depth_R,
                'p_raw':       p_raw,
                'low_depth':   low_depth,
            })

        # -- mutagenesis index -------------------------------------------------
        # total non-WT count = depth - WT count  (with pseudocount)
        wt_count_X   = int(x[wt_base])
        wt_count_R   = int(r[wt_base])
        nonwt_count_X = depth_X - wt_count_X
        nonwt_count_R = depth_R - wt_count_R

        # freq of any mutation in X and R (pseudocount on non-WT total)
        freq_mut_X = (nonwt_count_X + PSEUDOCOUNT) / (depth_X + PSEUDOCOUNT)
        freq_mut_R = (nonwt_count_R + PSEUDOCOUNT) / (depth_R + PSEUDOCOUNT)

        # WT-normalized mutagenesis index
        # = [freq_mut_X / freq_mut_R] x [freq_wt_R / freq_wt_X]
        mut_index     = (freq_mut_X / freq_mut_R) * wt_ratio
        log2_mut_index = np.log2(mut_index)

        # z-test: is total mutation frequency different between X and R?
        try:
            _, p_mi = proportions_ztest(
                count=[nonwt_count_X + PSEUDOCOUNT,
                       nonwt_count_R + PSEUDOCOUNT],
                nobs =[depth_X + PSEUDOCOUNT,
                       depth_R + PSEUDOCOUNT],
                alternative='two-sided'
            )
        except Exception:
            p_mi = np.nan

        mi_rows.append({
            'pos':             pos,
            'wt_base':         wt_base,
            'mut_index':       mut_index,
            'log2_mut_index':  log2_mut_index,
            'nonwt_count_X':   nonwt_count_X,
            'depth_X':         depth_X,
            'nonwt_count_R':   nonwt_count_R,
            'depth_R':         depth_R,
            'p_mi_raw':        p_mi,
            'low_depth':       low_depth,
        })

    return pd.DataFrame(fc_rows), pd.DataFrame(mi_rows)


# -- statistics ----------------------------------------------------------------

def compute_stats(rep_fc_list: list) -> tuple:
    """
    Merge FC and mutagenesis index across replicate pairs.
    Combine p-values with Stouffer's method, apply BH FDR correction.

    Returns:
        stats_df   — per-base FC statistics (one row per pos x mut_base)
        mi_df      — per-position mutagenesis index statistics
    """
    fc_frames = []
    mi_frames = []

    for rep_label, (fc_df, mi_df) in rep_fc_list:
        # -- per-base FC -------------------------------------------------------
        tmp = fc_df[['pos', 'wt_base', 'mut_base',
                     'log2fc', 'p_raw', 'low_depth',
                     'count_mut_X', 'depth_X',
                     'count_mut_R', 'depth_R']].copy()
        tmp = tmp.rename(columns={
            'log2fc':      f'log2fc_rep{rep_label}',
            'p_raw':       f'p_rep{rep_label}',
            'low_depth':   f'low_depth_{rep_label}',
            'count_mut_X': f'count_X_rep{rep_label}',
            'depth_X':     f'depth_X_rep{rep_label}',
            'count_mut_R': f'count_R_rep{rep_label}',
            'depth_R':     f'depth_R_rep{rep_label}',
        })
        fc_frames.append(tmp)

        # -- mutagenesis index -------------------------------------------------
        mf = mi_df[['pos', 'wt_base',
                    'log2_mut_index', 'p_mi_raw', 'low_depth',
                    'nonwt_count_X', 'depth_X',
                    'nonwt_count_R', 'depth_R']].copy()
        mf = mf.rename(columns={
            'log2_mut_index': f'log2_mi_rep{rep_label}',
            'p_mi_raw':       f'p_mi_rep{rep_label}',
            'low_depth':      f'low_depth_{rep_label}',
            'nonwt_count_X':  f'nonwt_X_rep{rep_label}',
            'depth_X':        f'depth_X_rep{rep_label}',
            'nonwt_count_R':  f'nonwt_R_rep{rep_label}',
            'depth_R':        f'depth_R_rep{rep_label}',
        })
        mi_frames.append(mf)

    # -- merge FC across reps --------------------------------------------------
    merged_fc = fc_frames[0]
    for f in fc_frames[1:]:
        merged_fc = merged_fc.merge(f, on=['pos', 'wt_base', 'mut_base'],
                                    how='inner')

    low_cols    = [c for c in merged_fc.columns if c.startswith('low_depth_')]
    log2fc_cols = [c for c in merged_fc.columns if c.startswith('log2fc_rep')]
    p_cols      = [c for c in merged_fc.columns if c.startswith('p_rep')]

    merged_fc['low_depth_flag'] = merged_fc[low_cols].any(axis=1)
    merged_fc = merged_fc.drop(columns=low_cols)
    merged_fc['mean_log2fc'] = merged_fc[log2fc_cols].mean(axis=1)
    merged_fc['p_combined']  = merged_fc.apply(
        lambda row: _stouffer(row, p_cols), axis=1)

    valid = ~merged_fc['p_combined'].isna()
    p_adj = np.full(len(merged_fc), np.nan)
    if valid.sum() > 0:
        _, p_adj[valid.values], _, _ = multipletests(
            merged_fc.loc[valid, 'p_combined'], method='fdr_bh')
    merged_fc['p_adj'] = p_adj

    fc_col_order = (['pos', 'wt_base', 'mut_base']
                    + log2fc_cols
                    + ['mean_log2fc', 'p_combined', 'p_adj', 'low_depth_flag']
                    + [c for c in merged_fc.columns
                       if c.startswith('count_') or c.startswith('depth_')])
    merged_fc = merged_fc[fc_col_order].sort_values(
        ['pos', 'mut_base']).reset_index(drop=True)

    # -- merge mutagenesis index across reps -----------------------------------
    merged_mi = mi_frames[0]
    for f in mi_frames[1:]:
        merged_mi = merged_mi.merge(f, on=['pos', 'wt_base'], how='inner')

    low_cols_mi  = [c for c in merged_mi.columns if c.startswith('low_depth_')]
    log2mi_cols  = [c for c in merged_mi.columns if c.startswith('log2_mi_rep')]
    p_mi_cols    = [c for c in merged_mi.columns if c.startswith('p_mi_rep')]

    merged_mi['low_depth_flag']  = merged_mi[low_cols_mi].any(axis=1)
    merged_mi = merged_mi.drop(columns=low_cols_mi)
    merged_mi['mean_log2_mi']    = merged_mi[log2mi_cols].mean(axis=1)
    merged_mi['p_mi_combined']   = merged_mi.apply(
        lambda row: _stouffer(row, p_mi_cols), axis=1)

    valid_mi = ~merged_mi['p_mi_combined'].isna()
    p_adj_mi = np.full(len(merged_mi), np.nan)
    if valid_mi.sum() > 0:
        _, p_adj_mi[valid_mi.values], _, _ = multipletests(
            merged_mi.loc[valid_mi, 'p_mi_combined'], method='fdr_bh')
    merged_mi['p_adj_mi'] = p_adj_mi

    mi_col_order = (['pos', 'wt_base']
                    + log2mi_cols
                    + ['mean_log2_mi', 'p_mi_combined', 'p_adj_mi',
                       'low_depth_flag']
                    + [c for c in merged_mi.columns
                       if c.startswith('nonwt_') or c.startswith('depth_')])
    merged_mi = merged_mi[mi_col_order].sort_values('pos').reset_index(drop=True)

    return merged_fc, merged_mi


def _stouffer(row, p_cols: list) -> float:
    """Combine p-values across replicates using Stouffer's method."""
    pvals = [row[c] for c in p_cols if not np.isnan(row[c])]
    if len(pvals) == 0:
        return np.nan
    if len(pvals) == 1:
        return pvals[0]
    _, p_combined = combine_pvalues(pvals, method='stouffer')
    return p_combined


# -- volcano plot --------------------------------------------------------------

def plot_volcano(stats_df: pd.DataFrame, title: str, out_path: str,
                 fc_thresh: float = 1.5, p_thresh: float = 0.05):
    df = stats_df.dropna(subset=['p_adj']).copy()
    df['neg_log10_p'] = -np.log10(df['p_adj'].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 5.5), dpi=150)

    fc_line = np.log2(fc_thresh)
    for xv in [-fc_line, fc_line]:
        ax.axvline(xv, color='#888', lw=0.8, ls='--', zorder=1)
    ax.axhline(-np.log10(p_thresh), color='#888', lw=0.8, ls='--', zorder=1)
    ax.axvline(0, color='#ccc', lw=0.6, zorder=1)

    for mut_base, grp in df.groupby('mut_base'):
        color = BASE_COLORS.get(mut_base, '#999')
        hi = grp[~grp['low_depth_flag']]
        lo = grp[ grp['low_depth_flag']]
        ax.scatter(hi['mean_log2fc'], hi['neg_log10_p'],
                   c=color, s=18, alpha=0.75, lw=0,
                   label=f'->{mut_base}', zorder=2)
        if len(lo):
            ax.scatter(lo['mean_log2fc'], lo['neg_log10_p'],
                       facecolors='none', edgecolors=color,
                       s=18, alpha=0.55, lw=0.8, zorder=2)

    sig = df[
        (df['p_adj'] < p_thresh) &
        (df['mean_log2fc'].abs() > fc_line) &
        (~df['low_depth_flag'])
    ].copy()
    sig['label'] = sig['wt_base'] + sig['pos'].astype(str) + sig['mut_base']
    for _, row in sig.nlargest(20, 'neg_log10_p').iterrows():
        ax.annotate(row['label'],
                    xy=(row['mean_log2fc'], row['neg_log10_p']),
                    xytext=(4, 2), textcoords='offset points',
                    fontsize=5.5, color='#222',
                    arrowprops=dict(arrowstyle='-', lw=0.4, color='#aaa'))

    ax.set_xlabel('mean log2(WT-normalized FC)', fontsize=11)
    ax.set_ylabel('-log10(p adj, BH)', fontsize=11)
    ax.set_title(title, fontsize=12, fontweight='bold')

    handles = [mpatches.Patch(facecolor=BASE_COLORS[b], label=f'->{b}')
               for b in BASES]
    handles.append(mlines.Line2D([0], [0], marker='o', color='w',
                                 markerfacecolor='none', markeredgecolor='#555',
                                 markersize=6, label='low depth'))
    ax.legend(handles=handles, title='mut base', fontsize=8,
              title_fontsize=8, frameon=False, loc='upper left')

    enr = sig[sig['mean_log2fc'] >  fc_line]
    dep = sig[sig['mean_log2fc'] < -fc_line]
    ax.text(0.98, 0.98, f'enriched: {len(enr)}',
            transform=ax.transAxes, ha='right', va='top',
            fontsize=8, color='#c0392b')
    ax.text(0.02, 0.98, f'depleted: {len(dep)}',
            transform=ax.transAxes, ha='left', va='top',
            fontsize=8, color='#2980b9')

    plt.tight_layout()
    plt.savefig(out_path, bbox_inches='tight',
                dpi=200 if out_path.endswith('.png') else None)
    print(f"Saved plot: {out_path}")
    plt.close()


# -- CLI -----------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--readcount-dir', required=True)
    p.add_argument('--virus',         required=True, choices=['S', 'B'])
    p.add_argument('--out-tsv',       default='wtnorm_fc.tsv')
    p.add_argument('--out-mi-tsv',    default='mutindex.tsv',
                   help='Output TSV for per-position mutagenesis index')
    p.add_argument('--out-plot',      default='wtnorm_volcano.pdf')
    p.add_argument('--title',         default='WT-normalized mutation FC')
    p.add_argument('--fc-thresh',     type=float, default=1.5)
    p.add_argument('--p-thresh',      type=float, default=0.05)
    return p.parse_args()


def main():
    args = parse_args()

    paired = find_sample_files(args.readcount_dir, args.virus)
    print(f"Found {len(paired)} replicate pair(s): reps {sorted(paired)}")

    rep_fc_list = []
    for rep in sorted(paired):
        pair  = paired[rep]
        x_df  = parse_readcount_file(pair['X'])
        r_df  = parse_readcount_file(pair['R'])
        fc_df, mi_df = wt_normalized_fc(x_df, r_df)
        rep_fc_list.append((rep, (fc_df, mi_df)))

    print("Computing statistics...")
    result_fc, result_mi = compute_stats(rep_fc_list)

    # save per-base FC TSV
    result_fc.to_csv(args.out_tsv, sep='\t', index=False, float_format='%.6f')
    n_sig = (result_fc['p_adj'] < args.p_thresh).sum()
    print(f"Saved: {args.out_tsv}  ({len(result_fc)} mutations, "
          f"{n_sig} significant at p_adj < {args.p_thresh})")

    # save mutagenesis index TSV
    result_mi.to_csv(args.out_mi_tsv, sep='\t', index=False,
                     float_format='%.6f')
    print(f"Saved: {args.out_mi_tsv}  ({len(result_mi)} positions)")

    print("Plotting volcano...")
    plot_volcano(result_fc, title=args.title, out_path=args.out_plot,
                 fc_thresh=args.fc_thresh, p_thresh=args.p_thresh)


if __name__ == '__main__':
    main()
