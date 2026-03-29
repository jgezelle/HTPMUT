#!/usr/bin/env python3
"""
wtnorm_fc.py  (reads directly from bam-readcount output)

WT-normalized fold change per mutation per position, across paired replicates.

Formula (per replicate pair, per non-WT base B at position i):
    FC(B,i) = [freq(B,X) / freq(B,R)] x [freq(WT,R) / freq(WT,X)]

Statistics:
    For each mutation, use a 2-proportion z-test on raw read counts
    (count of mut base vs depth) between X and R within each replicate.
    P-values are combined across replicates using Stouffer's method.
    Combined p-values are BH-corrected across all mutations.

    This is statistically appropriate because:
    - It uses raw counts directly, fully exploiting read depth
    - It tests whether mut frequency in X differs from R
    - It does not assume variance structure across replicates

Usage:
    python wtnorm_fc.py \\
        --readcount-dir  /path/to/05_readcount/260327 \\
        --virus          S \\
        --out-tsv        scnv_wtnorm_fc.tsv \\
        --out-plot       scnv_volcano.pdf \\
        --title          "scnv XR2"

    # --virus S  ->  SCNV samples
    # --virus B  ->  BENY samples
    # Replicates paired automatically: X1<->R1, X2<->R2

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
from scipy import stats
from scipy.stats import combine_pvalues
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportions_ztest

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'

EXCLUDE_START   = 29      # exclude positions 1-29 (leader / before stop site)
MIN_DEPTH_RATIO = 0.5     # minimum X/R depth ratio — catches transition zone
MIN_DEPTH_X   = 500   # require meaningful X coverage
MIN_DEPTH_R   = 1000  # R should always have deep coverage
BASES       = ['A', 'C', 'G', 'T']
PSEUDOCOUNT = 0.5
BASE_COLORS = {'A': '#E64B35', 'C': '#4DBBD5', 'G': '#F39B7F', 'T': '#3C5488'}


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
    # frequencies with pseudocount
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
        cond  = parts[3]    # X or R
        rep   = parts[4]    # 1 or 2
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
                     r_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute WT-normalized FC and 2-proportion z-test p-value for each
    non-WT substitution at each shared position.

    Returns long DataFrame:
        pos  wt_base  mut_base  fc  log2fc
        count_mut_X  depth_X  count_mut_R  depth_R
        p_value_raw  low_depth
    """
    shared = x_df.index.intersection(r_df.index)
    rows = []
    for pos in shared:
        if pos <= EXCLUDE_START:
            continue
            
        wt_base = x_df.loc[pos, 'ref']
        if wt_base not in BASES:
            continue
        x = x_df.loc[pos]
        r = r_df.loc[pos]

        depth_X = int(x['depth'])
        depth_R = int(r['depth'])
        if pos <= EXCLUDE_START:
            continue
        low_depth = (int(x['depth']) < MIN_DEPTH_X) or (int(r['depth']) < MIN_DEPTH_R)

        freq_wt_X = x[f'freq_{wt_base}']
        freq_wt_R = r[f'freq_{wt_base}']
        wt_ratio  = freq_wt_R / freq_wt_X   # WT(R) / WT(X) — normalization term

        for mut_base in BASES:
            if mut_base == wt_base:
                continue   # skip WT — FC is trivially 1, not informative

            count_X = int(x[mut_base])
            count_R = int(r[mut_base])

            # WT-normalized fold change
            fc = (x[f'freq_{mut_base}'] / r[f'freq_{mut_base}']) * wt_ratio
            log2fc = np.log2(fc)

            # 2-proportion z-test on raw counts: is freq(mut) different in X vs R?
            # counts = [mut in X, mut in R], nobs = [depth_X, depth_R]
            try:
                _, p_raw = proportions_ztest(
                    count   = [count_X + PSEUDOCOUNT, count_R + PSEUDOCOUNT],
                    nobs    = [depth_X + 4 * PSEUDOCOUNT,
                               depth_R + 4 * PSEUDOCOUNT],
                    alternative='two-sided'
                )
            except Exception:
                p_raw = np.nan

            rows.append({
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
    return pd.DataFrame(rows)


# -- statistics ----------------------------------------------------------------

def compute_stats(rep_fc_list: list) -> pd.DataFrame:
    """
    Merge FC and p-values across replicate pairs.
    Combine p-values using Stouffer's method.
    Apply BH FDR correction on combined p-values.
    """
    n_reps = len(rep_fc_list)

    # Merge all reps on pos/wt_base/mut_base
    merged = None
    for rep_label, df in rep_fc_list:
        tmp = df[['pos', 'wt_base', 'mut_base',
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
        merged = tmp if merged is None else merged.merge(
            tmp, on=['pos', 'wt_base', 'mut_base'], how='inner')

    low_cols    = [c for c in merged.columns if c.startswith('low_depth_')]
    log2fc_cols = [c for c in merged.columns if c.startswith('log2fc_rep')]
    p_cols      = [c for c in merged.columns if c.startswith('p_rep')]

    merged['low_depth_flag'] = merged[low_cols].any(axis=1)
    merged = merged.drop(columns=low_cols)

    # Mean log2FC across replicates
    merged['mean_log2fc'] = merged[log2fc_cols].mean(axis=1)

    # Combine p-values across replicates using Stouffer's method
    def stouffer_combine(row):
        pvals = [row[c] for c in p_cols if not np.isnan(row[c])]
        if len(pvals) == 0:
            return np.nan
        if len(pvals) == 1:
            return pvals[0]
        _, p_combined = combine_pvalues(pvals, method='stouffer')
        return p_combined

    merged['p_combined'] = merged.apply(stouffer_combine, axis=1)

    # BH FDR correction
    valid = ~merged['p_combined'].isna()
    p_adj = np.full(len(merged), np.nan)
    if valid.sum() > 0:
        _, p_adj[valid.values], _, _ = multipletests(
            merged.loc[valid, 'p_combined'], method='fdr_bh')
    merged['p_adj'] = p_adj

    col_order = (['pos', 'wt_base', 'mut_base']
                 + log2fc_cols
                 + ['mean_log2fc', 'p_combined', 'p_adj', 'low_depth_flag']
                 + [c for c in merged.columns
                    if c.startswith('count_') or c.startswith('depth_')])

    return merged[col_order].sort_values(['pos', 'mut_base']).reset_index(drop=True)


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

    # Label top significant hits
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
    p.add_argument('--readcount-dir', required=True,
                   help='Directory containing *.readcount.txt files')
    p.add_argument('--virus', required=True, choices=['S', 'B'],
                   help='Virus code: S = scnv, B = beny')
    p.add_argument('--out-tsv',   default='wtnorm_fc.tsv')
    p.add_argument('--out-plot',  default='wtnorm_volcano.pdf')
    p.add_argument('--title',     default='WT-normalized mutation FC')
    p.add_argument('--fc-thresh', type=float, default=1.5)
    p.add_argument('--p-thresh',  type=float, default=0.05)
    return p.parse_args()


def main():
    args = parse_args()

    paired = find_sample_files(args.readcount_dir, args.virus)
    print(f"Found {len(paired)} replicate pair(s): reps {sorted(paired)}")

    rep_fc_list = []
    for rep in sorted(paired):
        pair = paired[rep]
        x_df  = parse_readcount_file(pair['X'])
        r_df  = parse_readcount_file(pair['R'])
        fc_df = wt_normalized_fc(x_df, r_df)
        rep_fc_list.append((rep, fc_df))

    print("Computing statistics...")
    result = compute_stats(rep_fc_list)

    result.to_csv(args.out_tsv, sep='\t', index=False, float_format='%.6f')
    n_sig = (result['p_adj'] < args.p_thresh).sum()
    print(f"Saved: {args.out_tsv}  ({len(result)} mutations, "
          f"{n_sig} significant at p_adj < {args.p_thresh})")

    print("Plotting volcano...")
    plot_volcano(result, title=args.title, out_path=args.out_plot,
                 fc_thresh=args.fc_thresh, p_thresh=args.p_thresh)


if __name__ == '__main__':
    main()