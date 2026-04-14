#!/usr/bin/env python3
"""
extract_mutseqs_partial.py
Extract multi-mutant sequences from Xrn1-treated (NOLEAD) BAMs for
covariance model building and covariation analysis.

KEY DIFFERENCE FROM extract_mutseqs.py:
    Does NOT require reads to fully span [pos_min, pos_max].
    Any mapped read is included if it has >= min_muts mutations
    within the window. Positions within the window not covered
    by the read are filled with the reference base.

    This captures far more reads — especially useful when reads
    don't all start at exactly pos_min (e.g. after leader removal).

Pipeline:
  1. For each X replicate BAM, scan ALL mapped reads
  2. Count mutations falling within [pos_min, pos_max]
  3. Keep reads with >= min_muts mutations in window
  4. Reconstruct sequence in reference space (ref base for uncovered positions)
  5. Intersect: keep sequences seen in BOTH replicates
  6. Deduplicate, sort by mutation burden (descending)
  7. Random sample up to max_seqs
  8. Output FASTA + metadata TSV

Output files (one set per threshold):
  {prefix}.minN.fasta        — sequences for MUSCLE/MAFFT/mLocaRNA
  {prefix}.minN.meta.tsv     — per-sequence metadata
  {prefix}.minN.summary.txt  — run statistics

Usage:
    python extract_mutseqs_partial.py \\
        --bam-x1   path/to/X1.NOLEAD.sorted.bam \\
        --bam-x2   path/to/X2.NOLEAD.sorted.bam \\
        --ref      path/to/scnv_xr2_wt.fa \\
        --pos-min  30 \\
        --pos-max  80 \\
        --out      scnv_mutseqs_core \\
        --min-muts 2 3 \\
        --max-seqs 10000

Notes:
    - Mutations counted only within [pos_min, pos_max]
    - Uncovered positions within window filled with reference base
    - Intersection: exact sequence string must appear in both replicates
    - Random seed fixed for reproducibility (--seed to change)
"""

import argparse
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam

MIN_MAPQ = 20
MIN_BQ   = 20


# ── reference loading ─────────────────────────────────────────────────────────

def load_reference(ref_fa: str, contig: str,
                   pos_min: int, pos_max: int) -> str:
    """
    Load reference sequence for the window [pos_min, pos_max] (1-based inclusive).
    Returns uppercase DNA string of length (pos_max - pos_min + 1).
    """
    ref    = pysam.FastaFile(ref_fa)
    refseq = ref.fetch(contig, pos_min - 1, pos_max).upper().replace('U', 'T')
    ref.close()
    win_len = pos_max - pos_min + 1
    if len(refseq) != win_len:
        raise SystemExit(
            f"Reference window length mismatch: expected {win_len}, got {len(refseq)}")
    return refseq


# ── per-read sequence extraction ──────────────────────────────────────────────

def extract_reads(bam_path: str, ref_fa: str,
                  pos_min: int, pos_max: int,
                  min_muts: int) -> dict:
    """
    Scan one BAM file. For each mapped read with MAPQ >= MIN_MAPQ:
      - Count substitution mutations within [pos_min, pos_max]
      - Positions in window not covered by read → reference base
      - Include read if n_muts >= min_muts

    No spanning requirement — reads that partially overlap the window
    are included, with uncovered positions filled by reference.

    Returns dict: {sequence_string: n_muts}
    """
    bam    = pysam.AlignmentFile(bam_path, 'rb')
    contig = bam.references[0]
    refseq = load_reference(ref_fa, contig, pos_min, pos_max)

    results = {}
    seen = skipped_unmapped = skipped_mapq = skipped_no_overlap = 0

    for read in bam.fetch(contig):
        seen += 1

        if read.is_unmapped:
            skipped_unmapped += 1
            continue
        if read.mapping_quality < MIN_MAPQ:
            skipped_mapq += 1
            continue

        # check read overlaps window at all (0-based)
        if read.reference_end <= (pos_min - 1) or \
           read.reference_start >= pos_max:
            skipped_no_overlap += 1
            continue

        # build sequence array initialized to reference bases
        seq_arr = list(refseq)
        n_muts  = 0
        qseq    = read.query_sequence
        qqual   = read.query_qualities

        for qpos, rpos, rbase in read.get_aligned_pairs(with_seq=True):
            if qpos is None or rpos is None:
                continue                        # skip insertions/deletions
            pos1 = rpos + 1                     # convert to 1-based
            if pos1 < pos_min or pos1 > pos_max:
                continue                        # outside window
            if qqual is not None and qqual[qpos] < MIN_BQ:
                continue                        # low base quality
            alt      = qseq[qpos].upper().replace('U', 'T')
            ref_base = refseq[pos1 - pos_min]
            if alt in ('A', 'C', 'G', 'T') and alt != ref_base:
                seq_arr[pos1 - pos_min] = alt   # record mutation
                n_muts += 1

        if n_muts < min_muts:
            continue

        seq_str = ''.join(seq_arr)

        # keep entry with highest mutation count
        if seq_str not in results or results[seq_str] < n_muts:
            results[seq_str] = n_muts

    bam.close()

    print(f"    reads seen:            {seen}")
    print(f"    skipped unmapped:      {skipped_unmapped}")
    print(f"    skipped low mapq:      {skipped_mapq}")
    print(f"    skipped no overlap:    {skipped_no_overlap}")
    print(f"    passing reads:         {len(results)}")

    return results


# ── intersection + deduplication ──────────────────────────────────────────────

def intersect_replicates(rep1: dict, rep2: dict) -> list:
    """
    Keep only sequences present in both replicates.
    Returns list of (seq, n_muts) sorted by n_muts descending.
    """
    shared = set(rep1.keys()) & set(rep2.keys())
    result = []
    for seq in shared:
        n_muts = max(rep1[seq], rep2[seq])
        result.append((seq, n_muts))
    result.sort(key=lambda x: (-x[1], x[0]))
    return result


# ── mutation position annotation ──────────────────────────────────────────────

def get_mut_positions(seq: str, refseq: str, pos_min: int) -> str:
    """
    Return comma-separated list of mutation positions (1-based genomic).
    """
    muts = []
    for i, (s, r) in enumerate(zip(seq, refseq)):
        if s != r:
            muts.append(str(pos_min + i))
    return ','.join(muts)


# ── output writers ────────────────────────────────────────────────────────────

def write_outputs(sequences: list, refseq: str, pos_min: int,
                  prefix: str, min_muts: int,
                  max_seqs: int, seed: int):
    """
    Write FASTA and metadata TSV for one mutation threshold.
    """
    rng = random.Random(seed)
    if len(sequences) > max_seqs:
        sampled      = rng.sample(sequences, max_seqs)
        sampled.sort(key=lambda x: -x[1])
        n_sampled    = max_seqs
        was_sampled  = True
    else:
        sampled      = sequences
        n_sampled    = len(sequences)
        was_sampled  = False

    fasta_path = f"{prefix}.min{min_muts}.fasta"
    meta_path  = f"{prefix}.min{min_muts}.meta.tsv"
    summ_path  = f"{prefix}.min{min_muts}.summary.txt"

    # FASTA
    with open(fasta_path, 'w') as fh:
        for idx, (seq, n_muts) in enumerate(sampled, 1):
            mut_pos = get_mut_positions(seq, refseq, pos_min)
            fh.write(f">seq{idx:06d}_nmuts{n_muts}_{mut_pos}\n{seq}\n")
    print(f"  Wrote FASTA: {fasta_path}  ({n_sampled} sequences)")

    # metadata TSV
    with open(meta_path, 'w') as fh:
        fh.write("seq_id\tn_muts\tmut_positions\tsequence\n")
        for idx, (seq, n_muts) in enumerate(sampled, 1):
            mut_pos = get_mut_positions(seq, refseq, pos_min)
            seq_id  = f"seq{idx:06d}_nmuts{n_muts}"
            fh.write(f"{seq_id}\t{n_muts}\t{mut_pos}\t{seq}\n")
    print(f"  Wrote meta:  {meta_path}")

    # summary
    all_nmuts = [n for _, n in sequences]
    with open(summ_path, 'w') as fh:
        fh.write(f"min_muts_threshold:  {min_muts}\n")
        fh.write(f"sequences_in_intersection: {len(sequences)}\n")
        fh.write(f"sequences_output:    {n_sampled}\n")
        fh.write(f"randomly_sampled:    {was_sampled}\n")
        fh.write(f"random_seed:         {seed}\n")
        fh.write(f"partial_spanning:    True (no full-span requirement)\n")
        if all_nmuts:
            fh.write(f"n_muts_min:          {min(all_nmuts)}\n")
            fh.write(f"n_muts_max:          {max(all_nmuts)}\n")
            fh.write(f"n_muts_mean:         {sum(all_nmuts)/len(all_nmuts):.2f}\n")
        fh.write(f"\nmutation_burden_distribution:\n")
        burden = Counter(all_nmuts)
        for k in sorted(burden):
            fh.write(f"  {k} muts: {burden[k]} sequences\n")
    print(f"  Wrote summary: {summ_path}")

    return n_sampled


# ── CLI ───────────────────────────────────────────────────────────────────────

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--bam-x1',   required=True,
                   help='X replicate 1 NOLEAD sorted BAM')
    p.add_argument('--bam-x2',   required=True,
                   help='X replicate 2 NOLEAD sorted BAM')
    p.add_argument('--ref',      required=True,
                   help='Reference FASTA')
    p.add_argument('--pos-min',  type=int, required=True,
                   help='Window start (1-based)')
    p.add_argument('--pos-max',  type=int, required=True,
                   help='Window end (1-based)')
    p.add_argument('--out',      default='mutseqs',
                   help='Output prefix')
    p.add_argument('--min-muts', type=int, nargs='+', default=[2, 3],
                   help='Mutation threshold(s) (default: 2 3)')
    p.add_argument('--max-seqs', type=int, default=10000,
                   help='Max sequences per output FASTA (default: 10000)')
    p.add_argument('--seed',     type=int, default=42,
                   help='Random seed for sampling (default: 42)')
    return p.parse_args()


def main():
    args = parse_args()

    bam_tmp = pysam.AlignmentFile(args.bam_x1, 'rb')
    contig  = bam_tmp.references[0]
    bam_tmp.close()
    refseq  = load_reference(args.ref, contig, args.pos_min, args.pos_max)

    print(f"Reference contig:  {contig}")
    print(f"Window:            {args.pos_min}-{args.pos_max} "
          f"({args.pos_max - args.pos_min + 1} nt)")
    print(f"Thresholds:        {args.min_muts} mutations")
    print(f"Partial spanning:  True (uncovered positions = reference base)")
    print()

    min_thresh = min(args.min_muts)

    print(f"Scanning rep1 BAM (min_muts={min_thresh})...")
    rep1 = extract_reads(args.bam_x1, args.ref,
                         args.pos_min, args.pos_max, min_thresh)
    print()

    print(f"Scanning rep2 BAM (min_muts={min_thresh})...")
    rep2 = extract_reads(args.bam_x2, args.ref,
                         args.pos_min, args.pos_max, min_thresh)
    print()

    print("Intersecting replicates...")
    intersected = intersect_replicates(rep1, rep2)
    print(f"  Sequences in rep1:         {len(rep1)}")
    print(f"  Sequences in rep2:         {len(rep2)}")
    print(f"  Sequences in intersection: {len(intersected)}")
    print()

    for min_muts in sorted(args.min_muts):
        filtered = [(seq, n) for seq, n in intersected if n >= min_muts]
        print(f"Threshold >= {min_muts} mutations: {len(filtered)} sequences")
        if len(filtered) == 0:
            print(f"  WARNING: no sequences pass this threshold — skipping")
            continue
        write_outputs(filtered, refseq, args.pos_min,
                      args.out, min_muts, args.max_seqs, args.seed)
        print()

    print("Done.")


if __name__ == '__main__':
    main()
