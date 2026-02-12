#!/usr/bin/env python3
# readcount_to_subcounts.py
#
# Convert bam-readcount outputs (*.readcount.txt) into a tidy TSV of substitution “features”
# where each feature is pos_ref>alt (e.g., 45_A>G).
#
# This script:
#  - parses sample metadata from your dash-delimited sample names (ZR/ZX/DR/DX + replicate)
#  - extracts A/C/G/T base counts per position from bam-readcount lines
#  - writes one row per substitution (alt != ref)
#  - optionally restricts to your degenerate windows:
#        ZIKV: 31–96
#        DENV: 31–94

import argparse
import glob
import os
import re
from pathlib import Path
from typing import Dict, Tuple, Optional

BASES = ["A", "C", "G", "T"]


def parse_sample_metadata(sample_name: str) -> Tuple[str, str, str, Optional[str]]:
    """
    sample_name example:
      260120-r0399-ZX-01-A3466_R1
    returns: (virus, condition, replicate, code)
      virus: ZIKV/DENV
      condition: R/X
      replicate: 01/02/03
      code: ZR/ZX/DR/DX or None
    """
    parts = sample_name.split("-")

    code = None
    for p in parts:
        if p in {"ZR", "ZX", "DR", "DX"}:
            code = p
            break

    if code is None:
        return ("Unknown", "Unknown", "Unknown", None)

    virus = "ZIKV" if code[0] == "Z" else "DENV"
    condition = code[1]  # 'R' or 'X'

    replicate = "Unknown"
    idx = parts.index(code)
    if idx + 1 < len(parts):
        replicate = parts[idx + 1]  # e.g. "01"/"02"/"03"

    return (virus, condition, replicate, code)


def parse_base_counts(fields) -> Dict[str, int]:
    """
    bam-readcount line is tab-delimited:
      contig pos ref depth =:... A:count:... C:count:... G:count:... T:count:... N:count:...
    We extract the integer count (2nd colon-separated item) from tokens like 'A:4276735:...'
    """
    counts = {b: 0 for b in BASES}
    for tok in fields:
        # tokens begin like "A:123:..."
        if len(tok) >= 3 and tok[1] == ":" and tok[0] in BASES:
            parts = tok.split(":")
            if len(parts) >= 2:
                try:
                    counts[tok[0]] = int(parts[1])
                except ValueError:
                    pass
    return counts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--readcount_dir", required=True, help="Directory containing *.readcount.txt")
    ap.add_argument("--out_tsv", required=True, help="Output TSV path")
    ap.add_argument("--min_depth", type=int, default=1, help="Keep positions with depth >= this")
    ap.add_argument("--restrict", action="store_true",
                    help="Restrict to degenerate windows (ZIKV 31-96; DENV 31-94)")
    ap.add_argument("--zikv_deg", default="31-96", help="ZIKV degenerate window, e.g. 31-96")
    ap.add_argument("--denv_deg", default="31-94", help="DENV degenerate window, e.g. 31-94")
    args = ap.parse_args()

    def parse_window(s: str) -> Tuple[int, int]:
        m = re.match(r"^\s*(\d+)\s*-\s*(\d+)\s*$", s)
        if not m:
            raise SystemExit(f"Bad window format: {s} (expected like 31-96)")
        return (int(m.group(1)), int(m.group(2)))

    zikv_deg_start, zikv_deg_end = parse_window(args.zikv_deg)
    denv_deg_start, denv_deg_end = parse_window(args.denv_deg)

    files = sorted(glob.glob(os.path.join(args.readcount_dir, "*.readcount.txt")))
    if not files:
        raise SystemExit(f"No *.readcount.txt found in {args.readcount_dir}")

    out_dir = os.path.dirname(args.out_tsv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    with open(args.out_tsv, "w") as out:
        out.write("\t".join([
            "sample", "code", "virus", "condition", "replicate",
            "contig", "pos", "ref", "alt", "alt_count", "depth", "feature"
        ]) + "\n")

        for fp in files:
            sample = Path(fp).name.replace(".readcount.txt", "")
            virus, condition, rep, code = parse_sample_metadata(sample)

            with open(fp, "r") as f:
                for line in f:
                    if not line.strip():
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) < 5:
                        continue

                    contig = fields[0]
                    pos = int(fields[1])
                    ref = fields[2].upper().replace("U", "T")
                    depth = int(fields[3])

                    if depth < args.min_depth:
                        continue

                    if args.restrict:
                        if virus == "ZIKV" and not (zikv_deg_start <= pos <= zikv_deg_end):
                            continue
                        if virus == "DENV" and not (denv_deg_start <= pos <= denv_deg_end):
                            continue
                        # if Unknown, drop under restrict
                        if virus not in {"ZIKV", "DENV"}:
                            continue

                    bc = parse_base_counts(fields)

                    for alt in BASES:
                        if alt == ref:
                            continue
                        alt_count = bc.get(alt, 0)
                        feature = f"{pos}_{ref}>{alt}"
                        out.write("\t".join(map(str, [
                            sample, code or "NA", virus, condition, rep,
                            contig, pos, ref, alt, alt_count, depth, feature
                        ])) + "\n")


if __name__ == "__main__":
    main()