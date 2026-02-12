#!/usr/bin/env python3
# bp_score_from_l2fc.py
#
# Given substitution_l2fc.tsv (with per-substitution log2FC), compute base-pairing scores
# using a dot-bracket structure + its sequence and coordinate span on the reference.
#
# For each base pair (i,j) in the dot-bracket, mapped to reference positions:
#   maintain set = substitutions that keep canonical pairing with partner base
#   disrupt set  = substitutions that break pairing
#
# Score:
#   BPscore = mean(LFC_maintain) - mean(LFC_disrupt)
#
# Output: bp_scores.tsv with one row per base pair (restricted to degenerate window)

import argparse
import pandas as pd
from typing import List, Tuple, Dict

CANON = {("A","T"),("T","A"),("G","C"),("C","G"),("G","T"),("T","G")}  # treat U as T
BASES = ["A","C","G","T"]


def dotbracket_pairs(dot: str) -> List[Tuple[int,int]]:
    stack = []
    pairs = []
    for idx, ch in enumerate(dot, start=1):
        if ch == "(":
            stack.append(idx)
        elif ch == ")":
            if not stack:
                raise ValueError("Unbalanced dot-bracket: too many ')'")
            i = stack.pop()
            pairs.append((i, idx))
    if stack:
        raise ValueError("Unbalanced dot-bracket: too many '('")
    return pairs


def is_canonical(b1: str, b2: str) -> bool:
    b1 = b1.upper().replace("U","T")
    b2 = b2.upper().replace("U","T")
    return (b1, b2) in CANON


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--l2fc_tsv", required=True, help="substitution_l2fc.tsv")
    ap.add_argument("--out_tsv", required=True, help="Output bp_scores.tsv")
    ap.add_argument("--virus", required=True, choices=["ZIKV","DENV"])
    ap.add_argument("--dotbracket", required=True, help="Dot-bracket string for the span")
    ap.add_argument("--seq", required=True, help="Sequence string matching dotbracket span (RNA ok)")
    ap.add_argument("--start_pos", type=int, required=True,
                    help="Reference position of first character in dotbracket span (e.g., 24 or 25)")
    ap.add_argument("--deg_start", type=int, required=True, help="Degenerate window start (ref coords)")
    ap.add_argument("--deg_end", type=int, required=True, help="Degenerate window end (ref coords)")
    ap.add_argument("--min_subs_per_class", type=int, default=1,
                    help="Require at least this many maintaining AND disrupting substitutions per pair")
    args = ap.parse_args()

    dot = args.dotbracket.strip()
    seq = args.seq.strip().upper().replace("U","T")

    if len(dot) != len(seq):
        raise SystemExit(f"Dot-bracket length ({len(dot)}) != sequence length ({len(seq)}).")

    df = pd.read_csv(args.l2fc_tsv, sep="\t")
    df = df[df["virus"] == args.virus].copy()

    l2fc_map: Dict[str, float] = dict(zip(df["feature"], df["log2fc"]))

    pairs_local = dotbracket_pairs(dot)

    rows = []
    for i_local, j_local in pairs_local:
        i_pos = args.start_pos + i_local - 1
        j_pos = args.start_pos + j_local - 1

        # keep only pairs fully within degenerate region
        if not (args.deg_start <= i_pos <= args.deg_end and args.deg_start <= j_pos <= args.deg_end):
            continue

        ref_i = seq[i_local - 1]
        ref_j = seq[j_local - 1]

        maintain = []
        disrupt = []

        # substitutions at i
        for alt in BASES:
            if alt == ref_i:
                continue
            feat = f"{i_pos}_{ref_i}>{alt}"
            if feat in l2fc_map:
                (maintain if is_canonical(alt, ref_j) else disrupt).append(l2fc_map[feat])

        # substitutions at j
        for alt in BASES:
            if alt == ref_j:
                continue
            feat = f"{j_pos}_{ref_j}>{alt}"
            if feat in l2fc_map:
                (maintain if is_canonical(ref_i, alt) else disrupt).append(l2fc_map[feat])

        if len(maintain) < args.min_subs_per_class or len(disrupt) < args.min_subs_per_class:
            continue

        mean_m = sum(maintain) / len(maintain)
        mean_d = sum(disrupt) / len(disrupt)
        bpscore = mean_m - mean_d

        rows.append({
            "virus": args.virus,
            "i_pos": i_pos,
            "j_pos": j_pos,
            "ref_i": ref_i,
            "ref_j": ref_j,
            "n_maintain": len(maintain),
            "n_disrupt": len(disrupt),
            "mean_l2fc_maintain": mean_m,
            "mean_l2fc_disrupt": mean_d,
            "bp_score": bpscore
        })

    out = pd.DataFrame(rows).sort_values(["i_pos","j_pos"])
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote {len(out)} base-pair scores to {args.out_tsv}")


if __name__ == "__main__":
    main()