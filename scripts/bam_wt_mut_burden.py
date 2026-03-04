#!/usr/bin/env python3
import sys
from collections import Counter
import pysam

# args: bam ref_fa pos_min pos_max out_prefix
bam_path, ref_fa, pos_min, pos_max, out_prefix = sys.argv[1:]
pos_min = int(pos_min); pos_max = int(pos_max)

bam = pysam.AlignmentFile(bam_path, "rb")
ref = pysam.FastaFile(ref_fa)
contig = bam.references[0]

ref_seq = ref.fetch(contig, pos_min-1, pos_max).upper().replace("U","T")
win_len = pos_max - pos_min + 1
if len(ref_seq) != win_len:
    raise SystemExit("Reference window length mismatch")

min_mapq = 20
min_bq = 20

burden = Counter()
seen = used = unmapped = lowmapq = 0

for r in bam.fetch(contig):
    seen += 1
    if r.is_unmapped:
        unmapped += 1
        continue
    if r.mapping_quality < min_mapq:
        lowmapq += 1
        continue

    muts = 0
    qseq = r.query_sequence
    qqual = r.query_qualities

    for qpos, rpos, rbase in r.get_aligned_pairs(with_seq=True):
        if qpos is None or rpos is None:
            continue
        pos1 = rpos + 1
        if pos1 < pos_min or pos1 > pos_max:
            continue
        if qqual is not None and qqual[qpos] < min_bq:
            continue
        ref_base = ref_seq[pos1 - pos_min]
        alt = qseq[qpos].upper().replace("U","T")
        if alt in ("A","C","G","T") and alt != ref_base:
            muts += 1

    burden[muts] += 1
    used += 1

bam.close(); ref.close()

# burden table
with open(out_prefix + ".burden.tsv", "w") as out:
    out.write("n_mutations\tcount\n")
    for k in sorted(burden):
        out.write(f"{k}\t{burden[k]}\n")

# summary (1 line)
wt = burden.get(0, 0)
wt_frac = (wt / used) if used else 0.0

with open(out_prefix + ".summary.tsv", "w") as out:
    out.write("sample\tcontig\tpos_min\tpos_max\treads_seen\treads_used\twt_reads\twt_fraction\tskipped_unmapped\tskipped_lowmapq\n")
    out.write(f"{out_prefix.split('/')[-1]}\t{contig}\t{pos_min}\t{pos_max}\t{seen}\t{used}\t{wt}\t{wt_frac:.6f}\t{unmapped}\t{lowmapq}\n")