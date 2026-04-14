# HTP MUT 
## analysis of depleted mutant xrRNA sequences following exonuclease treatment 
### Mar 27, 2026 
### Beny and Soybean Cyst Nematode virus 
## Interactive job (`srun --pty`)
```bash
srun --pty -t 06:00:00 -c 4 --mem=16G \
--chdir=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT \
/bin/bash -l
```

### symlink 

```bash
export SRC=/groups/as6282_gp/data_bkup/jgg2144/20260325_Aviti_150x2_JGGMK_HTPMUT
ls $SRC/*.fastq.gz
```
Note – will use R1 only. <br>

```bash 
export ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
export RAW="$ROOT/work/01_raw"
cd "$RAW"

### (from `RAW` dir):
```bash
ln -s $SRC/A3474_R1.fastq.gz  260325-r0414-B-R-1-A3474_R1.fastq.gz
ln -s $SRC/A3471_R1.fastq.gz  260325-r0414-B-X-1-A3471_R1.fastq.gz

ln -s $SRC/A3470_R1.fastq.gz  260325-r0414-B-R-2-A3470_R1.fastq.gz
ln -s $SRC/A3454_R1.fastq.gz  260325-r0414-B-X-2-A3454_R1.fastq.gz

ln -s $SRC/A3449_R1.fastq.gz  260325-r0415-S-R-1-A3449_R1.fastq.gz
ln -s $SRC/A3408_R1.fastq.gz  260325-r0415-S-X-1-A3408_R1.fastq.gz

ln -s $SRC/A3405_R1.fastq.gz  260325-r0415-S-R-2-A3405_R1.fastq.gz
ln -s $SRC/A3477_R1.fastq.gz  260325-r0415-S-X-2-A3477_R1.fastq.gz

```

**check symlinks worked:**<br>
```bash
ls -lh
```

# 02_trim
```bash 
TRIM=$ROOT/work/02_trim/260327 # change for date!
LOGS=$TRIM/logs
mkdir -p "$TRIM" "$LOGS"
```
```bash
END3=AAAGAAACAACAACAACAAC
LEAD5=GGATTAATATAAT
REQUIRE_3P=true 
```
requiring that the `END3` sequence is present; in this case, true, because sequence length without exo treatment is 133 nt < 150 nt . 
<br>

## 3' end trim: 
```bash
cd "$TRIM"
```
```bash

for SAMPLE in "$RAW"/260325-*_R1.fastq.gz; do
  BASENAME=$(basename "$SAMPLE" .fastq.gz)

  if [ "$REQUIRE_3P" = true ]; then
    # Require the 3′ adaptor
    cutadapt -j 8 \
      -n 3 \
      -O 10 \
      -e 0.1 \
      -a "$END3" \
      -m 60 \
      --discard-untrimmed \
      -o "$TRIM/${BASENAME}.insert.fastq.gz" \
      "$SAMPLE" \
      > "$LOGS/${BASENAME}.cutadapt.3p.txt"
  else
    # Trim if present, but keep reads without it
    cutadapt -j 8 \
      -n 3 \
      -O 10 \
      -e 0.1 \
      -a "$END3" \
      -m 60 \
      -o "$TRIM/${BASENAME}.insert.fastq.gz" \
      "$SAMPLE" \
      > "$LOGS/${BASENAME}.cutadapt.3p.txt"
  fi
done
```
notes: <br>
 -- when `--discard-untrimmed` makes sense: 
adaptor must be present for the read to be interpretable ; adaptor presence is required to define molecule completeness
vs when it would be wrong: future libraries with longer sequence inserts may not capture the 3′ end <br>
-- `-a AAAGAAACAACAACAACAAC`
means:
“Find this sequence anywhere in the read (unanchored 3′ adaptor), and remove it and everything after it.” <br>
-- `-m 60` is used because the stop sites for this library are known and no sequences should be smaller than this. <br>
-- `-n 3` means that cutadapt can remove the adaptor up to three times per read, which is rare, but helps with concatemerized reads. <br>


## 5' end sort: 
```bash
for READS in "$TRIM"/*.insert.fastq.gz; do
  BASE=$(basename "$READS" .insert.fastq.gz)

  # split by leader
  cutadapt -j 8 \
    -b "$LEAD5" \
    -e 0.12 \
    --action=none \
    -o "$TRIM/${BASE}.HASLEAD.fastq.gz" \
    --untrimmed-output "$TRIM/${BASE}.NOLEAD.fastq.gz" \
    "$READS" \
    > "$LOGS/${BASE}.split5p.txt"

  # count - FASTQ files always have 4 lines per read... count lines and divide by 4 with wc -l

  HASLEAD_COUNT=$(( $(zcat "$TRIM/${BASE}.HASLEAD.fastq.gz" | wc -l) / 4 ))
  NOLEAD_COUNT=$(( $(zcat "$TRIM/${BASE}.NOLEAD.fastq.gz"  | wc -l) / 4 ))

  echo "$BASE HASLEAD: $HASLEAD_COUNT"
  echo "$BASE NOLEAD : $NOLEAD_COUNT"
done
```
Notes: <br>
--`-b "$LEAD5"` ; -b = front adaptor, asking whether reads start with the universal 5' end
<br>
--`-e 0.12` - allows 1-2 mismatches in the leader sequence (OK for now, avoids false negatives) <br>
--`--action=none` = if the leader is found, it should not be trimmed, but it will be accounted for <br>
-- `-o "$TRIM/${BASE}.HASLEAD.fastq.gz"` ; these reads matched the leader, still contain the leader, expected for R samples, but for X samples indicative of incomplete Xrn1 loading <br>
-- `--untrimmed-output "$TRIM/${BASE}.NOLEAD.fastq.gz"` These reads did NOT match the leader, expected for X samples, not for R samples. This flag creates the split.<br>


### make a summary table of the leader splits 
```bash 
SCRIPTS="$ROOT/scripts"
nano "$SCRIPTS/leader_summary.sh" # change trim directory at top of script! 
bash "$SCRIPTS/leader_summary.sh"
```
Loops over every *.HASLEAD.fastq.gz file in **TRIM** ; Finds the matching .NOLEAD.fastq.gz for each sample; Counts reads (lines / 4); Computes total reads per sample; Computes fraction of reads retaining the leader (HASLEAD / Total); Writes all of this to leader_summary.tsv


# 03_qc
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM=$ROOT/work/02_trim/260327 # change for date
QC=$ROOT/work/03_qc/260327 # change for date
mkdir -p "$QC"

# FastQC on all trimmed reads 
fastqc -t 8 -o "$QC" "$TRIM"/*.insert.fastq.gz

# Optional: leader-split QC 
fastqc -t 8 -o "$QC" "$TRIM"/*.HASLEAD.fastq.gz
fastqc -t 8 -o "$QC" "$TRIM"/*.NOLEAD.fastq.gz
```
```bash
# Aggregate into MultiQC report
multiqc -o "$QC" "$QC"
```
# 04_align
```bash 
cd "$ROOT/data/refseqs"
nano refseqs/scnv_xr2_wt.fa
nano refseqs/beny_xr_wt.fa
```
for now, both references include universal 5' end, but not the universal 3' end. <br>
**`scnv_xr2_wt.fa`** <br>
`GGATTAATATAATGGCGTGATTATTTAGCCCGTCAGCTTGACGTTAACTGCCACTTTGGTTGAAGTGTGATCAACCGTGCCTGGGGCGAGCATATCGGCCCAT` <br>
**`beny_xr_wt.fa`**<br>
`GGATTAATATAATATTTTATTTTCTTTTGGTGTAATCGTCCGAAGACGTTAAACTACACGTGATTTCACGGTGTTCGGTGAGAAGATTGTTTAACGGTGTTAC`


– notes - <br>
-will use `--very-sensitive-local` so that Bowtie2 can align part of the read, even if it isn't a perfect match (contains muts.)<br>
-recall unmatched ends will be soft-clipped <br>
–recall that ligation oligo added max 8 bases, which will be soft-clipped in either case <br>

**index** - create the binary lookup tables Bowtie2 uses for fast matching: <br>
```bash
cd "$ROOT/data/refseqs"
bowtie2-build scnv_xr2_wt.fa scnv_xr2_wt
bowtie2-build beny_xr_wt.fa beny_xr_wt
```
now, navigate to align dir: 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
mkdir -p "$ALIGN"
cd "$ALIGN"
```

### try align with one sample pair first: 
```bash
SAMPLE=260325-r0415-S-X-2-A3477_R1
REF="$ROOT/data/refseqs/scnv_xr2_wt"
```

```bash 
bowtie2 -p 8 \
  --very-sensitive-local \
  -x "$REF" \
  -U "$TRIM/${SAMPLE}.NOLEAD.fastq.gz" \
  2> "${SAMPLE}.bowtie2.log" \
| samtools view -bS - \
| samtools sort -o "${SAMPLE}.sorted.bam"

samtools index "${SAMPLE}.sorted.bam"
```
inspect (`flagstat` gives alignment information like mapped vs. unmapped reads):

```bash
samtools flagstat ${SAMPLE}.sorted.bam
less ${SAMPLE}.bowtie2.log
```

## loop for leader-split reads 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"
mkdir -p "$ALIGN"
cd "$ALIGN"
```
alignment loop for both leader-containing and no-lead samples after the leader-split step: <br>
``` bash
for READS in "$TRIM"/260325-*-R-*.HASLEAD.fastq.gz "$TRIM"/260325-*-X-*.NOLEAD.fastq.gz; do
  [[ -e "$READS" ]] || continue

  SAMPLE=$(basename "$READS" .fastq.gz)

  # Choose reference by virus label in filename
  if [[ "$SAMPLE" == *-S-* ]]; then
    REF="$REFDIR/scnv_xr2_wt"
  elif [[ "$SAMPLE" == *-B-* ]]; then
    REF="$REFDIR/beny_xr_wt"
  else
    echo "Unknown virus for $SAMPLE — skipping"
    continue
  fi

  echo "Aligning $SAMPLE"

  bowtie2 \
    -p 8 \
    --very-sensitive-local \
    -x "$REF" \
    -U "$READS" \
    2> "${SAMPLE}.bowtie2.log" \
  | samtools view -bS - \
  | samtools sort -o "${SAMPLE}.sorted.bam"

  samtools index "${SAMPLE}.sorted.bam"

  samtools flagstat "${SAMPLE}.sorted.bam" \
    > "${SAMPLE}.flagstat.txt"
done
```
```bash
 grep "overall alignment rate" *.bowtie2.log
```

 # 05_readcount 
 ## (and fc calculation and plotting)
 *can reorg later* <br>

```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"
SCRIPTS="$ROOT/scripts"
# run bam-readcount 
READCOUNT_DIR=$ROOT/work/05_readcount/260327
mkdir -p "$READCOUNT_DIR"


for BAM in "$ALIGN"/*.sorted.bam; do
  BASE=$(basename "$BAM" .sorted.bam)
  [[ "$BASE" == *-R-* && "$BASE" != *.HASLEAD ]] && continue
  [[ "$BASE" == *-X-* && "$BASE" != *.NOLEAD  ]] && continue
  [[ "$BASE" == *-S-* ]] && REF="$REFDIR/scnv_xr2_wt.fa"
  [[ "$BASE" == *-B-* ]] && REF="$REFDIR/beny_xr_wt.fa"
  IDXLINE=$(samtools idxstats "$BAM" | head -n 1)
  CONTIG=$(echo "$IDXLINE" | cut -f1)
  LEN=$(echo "$IDXLINE" | cut -f2)
  bam-readcount -f "$REF" -d 11000000 -q 20 -b 20 \
    "$BAM" "${CONTIG}:1-${LEN}" > "$READCOUNT_DIR/${BASE}.readcount.txt"
done
```
loops over every sorted BAM and for each one:

Filters by condition — skips R samples that aren't HASLEAD, and X samples that aren't NOLEAD <br>
Sets the reference — picks scnv or beny FASTA based on the S/B in the filename <br>
Auto-detects the contig region <br>
Runs bam-readcount — walks every position in the contig and for each one reports the exact count of every base (A, C, G, T) seen in reads passing the quality filter (MAPQ≥20, baseQ≥20). Output is one line per position. <br>

outputs a text file where every row = one genomic position, with exact base counts at that position. This then feeds the fold-change calculation <br>

### (Originally written to feed into DESeq2): 
(optional, for archiving) <br>
```bash
# build count matrices 
export READCOUNT_DIR && export OUTDIR=$ROOT/work/06_counts/260327
mkdir -p "$OUTDIR" && python "$SCRIPTS/readcount_to_matrix_v2.py"
```
 modified prev. version of `readcount_to_matrix.py` to accomodate new samples: 
 ```
 def parse_sample(sample):
    # e.g. 260325-r0415-S-X-2-A3477_R1.NOLEAD
    parts = sample.split("-")
    virus_code = parts[2]      # S or B
    condition  = parts[3]      # X or R
    replicate  = parts[4]      # 1 or 2

    virus = "SCNV" if virus_code == "S" else "BENY"
    return virus, condition, replicate
```

### FC calc (WT normalized): 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"
SCRIPTS="$ROOT/scripts"
READCOUNT_DIR=$ROOT/work/05_readcount/260327
PLOTS="$READCOUNT_DIR/plots"
mkdir -p "$PLOTS"

# WT-normalized FC + volcano — run once per virus

python "$SCRIPTS/wtnorm_fc.py" --readcount-dir $READCOUNT_DIR --virus S \
  --title "scnv XR2" --out-tsv "$PLOTS/scnv_wtnorm_fc.tsv" --out-plot "$PLOTS/scnv_volcano.pdf"

python "$SCRIPTS/wtnorm_fc.py" --readcount-dir $READCOUNT_DIR --virus B \
  --title "beny XR" --out-tsv "$PLOTS/beny_wtnorm_fc.tsv" --out-plot "$PLOTS/beny_volcano.pdf"
```
### Notes: 
### `wtnorm_fc.py`: <br> 
- reads .readcount.txt files directly (the **bam-readcount** file), auto-detects replicate pairs by parsing the filename (rep field pos. 4), and routes X/R pairing automatically 

calculation:
```
For each non-WT base B at position i, per replicate pair:

  freq(B, X)  = count(B in NOLEAD) / total depth at pos i in X
  freq(B, R)  = count(B in HASLEAD) / total depth at pos i in R
  freq(WT, X) = count(WT in NOLEAD) / total depth at pos i in X
  freq(WT, R) = count(WT in HASLEAD) / total depth at pos i in R

  FC = [freq(B,X) / freq(B,R)] × [freq(WT,R) / freq(WT,X)]
```
WT is considered at the normalization step — the second term [freq(WT,R) / freq(WT,X)] is what makes WT cancel to 1, because when B=WT those two terms are identical and FC=1 exactly. For any other base, asks "how much did this mutation change relative to how much WT itself changed at this position."

make sure the beginning "leader" + WT portion upstream of stop site is not included in the analysis. there is no expected X..NOLEAD coverage at those bases because they should be degraded. this exclusion is in the python script.

Example, if these bases are included: 
```
At position 28, X has ~466K reads and R has ~7.5M. Almost no reads in X means the counts for every non-WT base in X are near zero — just 0 or 1 reads. With pseudocount=1:

freq(mut, X) = (0 + 1) / (466429 + 4) ≈ 0.0000021
freq(mut, R) = (0 + 1) / (7474700 + 4) ≈ 0.00000013
So freq(mut, X) / freq(mut, R) ≈ 16
```
solution was to exclude those bases, as well as a depth ratio filter to directly address the cause <br> 
"don't trust this position because X and R have fundamentally different coverage, so the FC estimate is unreliable regardless of what the counts say" <br>

```
raw text file          →   list of dicts      →   DataFrame      →   DataFrame with frequencies
(one line per pos)         (one dict per pos)     (table)            (ready for FC calculation)
     ↑                           ↑                    ↑                        ↑
  open(path)               rows.append()        pd.DataFrame(rows)      freq_A = (A+pseudo)/depth
  
```bash
#use awk to check what the depth ratio looks like at the beginning of the sequence going into the stop site: 
awk 'NR==FNR{x[$2]=$4; next} ($2 in x){print $2, x[$2], $4, x[$2]/$4}' \
  $READCOUNT_DIR/260325-r0415-S-X-1-A3408_R1.NOLEAD.readcount.txt \
  $READCOUNT_DIR/260325-r0415-S-R-1-A3449_R1.HASLEAD.readcount.txt | head -45

# see that pos. 29-30 marks the change in depth threshold 
```
```
1 16 7176535 2.22949e-06
2 34 7350890 4.62529e-06
3 509 7413779 6.86559e-05
4 1204 7445885 0.0001617
5 2658 7459136 0.000356342
6 3806 7454836 0.000510541
7 5183 7453716 0.000695358
8 8016 7471754 0.00107284
9 11172 7425112 0.00150462
10 15419 7464611 0.00206561
11 17418 7452170 0.00233731
12 18871 7450133 0.00253297
13 22046 7475407 0.00294914
14 22886 7463070 0.00306657
15 23852 7408890 0.00321938
16 28472 7468008 0.00381253
17 28614 7343993 0.00389625
18 29677 7476124 0.00396957
19 29286 7014862 0.00417485
20 32858 7446641 0.00441246
21 36492 7464602 0.00488867
22 45423 7408921 0.00613085
23 58859 7467103 0.00788244
24 83380 7466818 0.0111667
25 109260 7466020 0.0146343
26 134639 7461952 0.0180434
27 267513 7465072 0.0358353
28 466429 7474700 0.062401
29 3278126 7476026 0.438485
30 6212880 7410176 0.838425
31 7842040 7456386 1.05172
32 7853107 7419620 1.05842
33 7929160 7463285 1.06242
34 7910143 7474602 1.05827
35 7923792 7453887 1.06304
36 7939135 7454756 1.06498
37 7919397 7450673 1.06291
38 7934140 7465211 1.06282
39 7934549 7447008 1.06547
40 7942639 7465681 1.06389
41 7919225 7448956 1.06313
42 7980715 7478661 1.06713
43 7905138 7407521 1.06718
44 7941841 7445941 1.0666
45 7957026 7445320 1.06873
```

still might be worth excluding start 30 which *just* makes the coverage ratio threshold. Could change threshold to 0.9 (pos 30 @ 0.838425)

### improvements to plot 
`plot_volcano.py` 

```bash
# pip install adjustText

python "$SCRIPTS/plot_volcano.py" \
  --tsv "$PLOTS/scnv_wtnorm_fc.tsv" \
  --title "scnv XR2" \
  --out "$PLOTS/scnv"
```
Notes - this needs a lot of work, very cluttered, don't use for now. Also probably best to do volcano plot in R. <br>

### heat map plot 
### `plot_heatmap.py`

```bash
python "$SCRIPTS/plot_heatmap.py" \
  --tsv   "$PLOTS/scnv_wtnorm_fc.tsv" \
  --title "scnv XR2" \
  --out   "$PLOTS/scnv_heatmap.pdf"
```
Notes on design:
- Non-significant cells are faded (18% opacity) so the significant signal stands out clearly (pass `--no-mask` to show everything at full color)
- WT base outlined in black per position from the wt_base column directly
`--vmax 4.0` by default (most depleted hits are around −7, so they'll be saturated blue, ... care about the pattern not the exact magnitude. Lower to `--vmax 2.0` for gradient in the milder hits
`--region 30 60` - zoom into the core structure without rerunning anything
`--structure` is fully optional — omit it and the structure panel just doesn't appear

updates needed: 
- fix overlapping text 
- make room along x axis for position (base #), and WT NT identity at that position, so the WT sequence aligns along the heat map squares and depleted mutations can be directly evaluated 
- ensure statistical method makes sense (p value means... )
  - What the p-value means here: the 2-proportion z-test asks "is the frequency of base B at position i significantly different between X and R?" The Stouffer combination asks "is this difference consistent across replicates?" The BH correction controls false discovery rate across all tested mutations. So p_adj < 0.05 means: this mutation's frequency difference between treated and control is unlikely to be due to sampling noise, after correcting for multiple testing. It does NOT account for the WT normalization step statistically — the WT normalization is applied to the FC only, the p-value is purely on the raw proportion difference
  - How replicates and depth are considered: replicates are handled by computing independent z-tests per replicate pair and combining with Stouffer's method.
  - ?? "Read depth is implicitly handled by the z-test — at 7M reads, even tiny frequency differences are significant, which is why so many hits saturate." ISSUE - statistical significance is essentially guaranteed at this depth, so the biologically meaningful axis is really the effect size (log2FC), not the p-value. The significance mask in the heatmap may not be adding much information. (what kind of statistical test WOULD BE BIOLOGICAL?)? 
- how are replicates and read depth considered
- parallel approach in R? 
- how else input can be used to make plots or even go into DE analysis ? 

# edits to plot script (v2)

Key changes:
- WT base row is now a proper colored tile row below the heatmap — same width as heatmap columns, base letter inside each tile (A=red, C=cyan, G=green, T=navy), position number shown every 5th position inside the tile. 
- No significance mask — full color everywhere, --vmax controls the saturation point. Start with --vmax 6 given the strongest hits are around −7, which will show gradient in the milder hits while saturating the core depleted positions.
Colorbar label now explicitly states what the statistic is and what direction means what.
- `--vmax` is now clearly documented — values beyond ±vmax are shown at the saturated color, not clipped to NaN.

```bash
python "$SCRIPTS/plot_heatmap_v2.py" \
  --tsv   "$PLOTS/scnv_wtnorm_fc.tsv" \
  --title "scnv XR2" \
  --vmax  6 \
  --out   "$PLOTS/scnv_heatmap.pdf"
```

Also made an R version, can run after scp to local: 
```bash
Rscript plot_heatmap_v2.R \
  --tsv scnv_wtnorm_fc.tsv \
  --title "scnv XR2" \
  --vmax 6 \
  --out scnv_heatmap.pdf
```

### aesthetic changes (bigger font and better base colors)

- A = Teal (#008B8B), 
- C = Maroon (#800000), 
- G = Emerald (#2E8B57), 
- T = Gold (#B8860B)

```bash
python "$SCRIPTS/plot_heatmap_v3.py" \
  --tsv   "$PLOTS/scnv_wtnorm_fc.tsv" \
  --title "scnv XR2" \
  --vmax  6 \
  --out   "$PLOTS/scnv_heatmap.pdf"
```
**summary on visualization vs. stats so far**: <br>
What the heatmap shows:

`mean_log2fc` = simple average of  `log2fc_rep1` and `log2fc_rep2`<br>
(no p-value, no weighting by depth, solely mean effect size)

What the `wtnorm_fc.py` script computed (in the TSV):<br>
- Per replicate: independent **2-proportion z-test** (X counts vs R counts at each position)
- Across replicates: Stouffer combination of the two p-values – `p_combined`
- Multiple testing correction – `p_adj`
- Both replicates contribute equally to `mean_log2fc`

**Ways the replicates are reflected:** <br>

- The mean log2FC averages both reps, so a position that's strongly depleted in both reps will show strong blue; if only one rep shows it, it gets diluted <br>
- The p-values in the TSV capture replicate consistency — but those are not currently visualized in the heatmap. <br>

- the heatmap is purely an effect size visualization.
- TSV shows the statistics <br>
- consider adding a dot in the p_adj < 0.05 cells <br>

- potential statistical issues: 
  - Z score is the frequency difference in X vs. R over the standard error 
  - at 7 million reads, the standard error on any frequency estimate is tiny — around 1/√7,000,000 ≈ 0.00038
```
freq difference = 0.0002
SE              ≈ 0.00038 / √2 ≈ 0.00027 (2 samples assumption.. )
z-score         = 0.0002 / 0.00027 ≈ 0.74   → p ≈ 0.46   (not significant)
freq difference = 0.001
z-score         ≈ 3.7    → p ≈ 0.0002 (significant)
```


for now, keep BH correction, use it to exclude genuinely noisy low-count mutations, but interpret the heatmap by effect size <br>

notes - 
- try to combine plot with shapemapper DMS reactivity? <br>
walk through code for `wtnorm_fc` and `plot_heatmap_v3` scripts. 
- additional output, PSL (Pattern Space Layout) file to store pairwise sequence alignments? can use `sam2psl`
- next / supplemental data: just the one hit reads, or sequence-level analysis 
- optimize normalization (cut off 5' end?)

# edits on Mar 31, 2026 
## mutagenesis index that can be superimposed on secondary structure 

A "mutagenesis index" — the aggregate mutation tolerance at each position regardless of which base it mutates to —
- collapses 3 mutations into one number per position. 
- Positions that are structurally critical will be depleted for ALL substitutions

The math: <br>
For each position i, per replicate:
```
mut_index(i) = sum of freq(all non-WT bases, X) / sum of freq(all non-WT bases, R)
               × [freq(WT, R) / freq(WT, X)]
```
Which simplifies to:
```
mut_index(i) = [1 - freq(WT, X)] / [1 - freq(WT, R)]  ×  [freq(WT, R) / freq(WT, X)]
```

total mutation frequency in treated vs control, WT-normalized. 
A value < 1 (log2 < 0) means mutations are depleted at that position — it's structurally intolerant. A value ≈ 1 means neutral.

So, **one number** per position rather than 3, which maps directly onto a secondary structure diagram as a color or thickness on each nucleotide.

### `plot_mutindex.py`
### `wtnorm_fc_aggreg`
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"
SCRIPTS="$ROOT/scripts"
READCOUNT_DIR=$ROOT/work/05_readcount/260327
PLOTS="$READCOUNT_DIR/plots"

# rerun wtnorm_fc.py, updated (wtnorm_fc_aggreg) — now outputs both TSVs
python "$SCRIPTS/wtnorm_fc_aggreg.py" \
  --readcount-dir $READCOUNT_DIR \
  --virus S \
  --out-tsv    "$PLOTS/scnv_wtnorm_fc_aggreg.tsv" \
  --out-mi-tsv "$PLOTS/scnv_mutindex.tsv" \
  --out-plot   "$PLOTS/scnv_volcano_aggreg.pdf" \
  --title "scnv XR2"

# 2. plot mutagenesis index
python "$SCRIPTS/plot_mutindex.py" \
  --tsv   "$PLOTS/scnv_mutindex.tsv" \
  --title "scnv XR2" \
  --out   "$PLOTS/scnv_mutindex" \
  --vmax  6
```

This produces three files: `scnv_mutindex_barplot.pdf`, `scnv_mutindex_heatmap.pdf`, and `scnv_mutindex.varna`.

**The VARNA file** contains one line per position with the log2 mutagenesis index value. 

*In VARNA: open structure, then go to **Annotations → Color map → Load** and select the `.varna` file. Set the color scale range to match `--vmax` (e.g. −6 to +6) with a blue-white-red diverging palette — blue nucleotides are structurally intolerant positions, white are neutral.*

**The math for the mutagenesis index** collapses all three non-WT bases into one number per position:
```
MI(i) = [freq(any mutation, X) / freq(any mutation, R)] × [freq(WT,R) / freq(WT,X)]
       = [(1 - freq_WT_X) / (1 - freq_WT_R)] × [freq_WT_R / freq_WT_X]
```


# Covariation downstream analysis on Apr 3, 2026 
## 06_mutseqs
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT

TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"
SCRIPTS="$ROOT/scripts"
READCOUNT_DIR=$ROOT/work/05_readcount/260327

ALIGN="$ROOT/work/04_align/260327"
MUTSEQ="$ROOT/work/06_mutseqs/260403"
mkdir -p "$MUTSEQ"

python "$SCRIPTS/extract_mutseqs.py" \
  --bam-x1  "$ALIGN/260325-r0415-S-X-1-A3408_R1.NOLEAD.sorted.bam" \
  --bam-x2  "$ALIGN/260325-r0415-S-X-2-A3477_R1.NOLEAD.sorted.bam" \
  --ref     "$REFDIR/scnv_xr2_wt.fa" \
  --pos-min 30 \
  --pos-max 101 \
  --out     "$MUTSEQ/scnv_mutseqs" \
  --min-muts 3 4 \
  --max-seqs 10000
```

- outputs six files — two sets (min3, min4) of FASTA + metadata TSV + summary. Key design decisions:

- Full-spanning requirement — `read.reference_start > (pos_min-1)` and `read.reference_end < pos_max` — only reads that cover the entire 30 to 101 window are included. - every sequence in the FASTA has complete information at every position.

- Intersection logic — exact sequence string match between rep1 and rep2. This is strict but appropriate — a sequence seen in both independent replicates is very likely a genuinely tolerated mutant rather than a PCR or sequencing artifact.

- FASTA header format — >seq000001_nmuts4_32,35,41,67 — encodes the sequence rank, mutation count, and exact positions mutated. This makes it easy to filter or annotate after alignment.

output: <br>

```
Scanning rep1 BAM (min_muts=3)...
    reads seen:          8057289
    skipped unmapped:    0
    skipped low mapq:    1822
    skipped non-spanning:1779810
    passing reads:       228381

Scanning rep2 BAM (min_muts=3)...
    reads seen:          5550074
    skipped unmapped:    0
    skipped low mapq:    2282
    skipped non-spanning:1371178
    passing reads:       238687

Intersecting replicates...
  Sequences in rep1:         228381
  Sequences in rep2:         238687
  Sequences in intersection: 83215

Threshold >= 3 mutations: 83215 sequences
  Wrote FASTA: /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min3.fasta  (10000 sequences)
  Wrote meta:  /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min3.meta.tsv
  Wrote summary: /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min3.summary.txt

Threshold >= 4 mutations: 2476 sequences
  Wrote FASTA: /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min4.fasta  (2476 sequences)
  Wrote meta:  /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min4.meta.tsv
  Wrote summary: /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs.min4.summary.txt
  ```

check the summary files: <br>
```bash
cat "$MUTSEQ/scnv_mutseqs.min3.summary.txt"
cat "$MUTSEQ/scnv_mutseqs.min4.summary.txt"
```
output: <br>
```
min_muts_threshold:  3
sequences_in_intersection: 83215
sequences_output:    10000
randomly_sampled:    True
random_seed:         42
n_muts_min:          3
n_muts_max:          5
n_muts_mean:         3.03

mutation_burden_distribution:
  3 muts: 80739 sequences
  4 muts: 2471 sequences
  5 muts: 5 sequences
min_muts_threshold:  4
sequences_in_intersection: 2476
sequences_output:    2476
randomly_sampled:    False
random_seed:         42
n_muts_min:          4
n_muts_max:          5
n_muts_mean:         4.00

mutation_burden_distribution:
  4 muts: 2471 sequences
  5 muts: 5 sequences
```
83k sequences with ≥3 mutations — the vast majority (80,739 / 97%) have exactly 3 mutations, with a sharp drop to 2,471 at 4 and only 5 at 5. This is expected from a degenerate library where mutation frequency is tuned low — most sequences that survive have the minimum number of mutations needed to be informative.
The min4 set (2,476 sequences) is a much better seed alignment for several reasons:

Small enough to align properly and inspect manually
4 mutations gives R-scape more covariation signal to detect real base pairs
The 3-mutation set at 83k would be dominated by sequences that are nearly identical to WT, diluting covariation signal
10k random sample from the 3-mutation set would be fine computationally but statistically weaker

The 5 sequences with 5 mutations are interesting — those are the most mutated tolerated sequences in the pool. Worth looking at them specifically:
```bash 
grep "nmuts5" "$MUTSEQ/scnv_mutseqs.min4.fasta"
```

## 07_alignment 
### use the minimum 4 mutants file first. 
```bash 
cd "$MUTSEQ"
ALIGNDIR="$MUTSEQ/07_align"
mkdir -p "$ALIGNDIR"

# MAFFT — fast, handles this size easily
mafft --auto \
  scnv_mutseqs.min4.fasta > "$ALIGNDIR/scnv_mutseqs.min4.afa"

# check alignment
grep -c "^>" scnv_mutseqs.min4.afa   # should be 2476
head -4 scnv_mutseqs.min4.afa
```
### instead of MAFFT, try direct `mLocaRNA` 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
SCRIPTS="$ROOT/scripts"
MUTSEQ="$ROOT/work/06_mutseqs/260403"
ALIGNDIR="$MUTSEQ/07_align"
mkdir -p "$ALIGNDIR"
cd "$ALIGNDIR"
```
check version 2 of locarna documentation: 
```bash
mlocarna --man 2>&1 > "$MUTSEQ/locarna_man.txt"
cat "$MUTSEQ/locarna_man.txt"
```
sequences are 72nt, 2,476 of them, substitutions only, so the "large instances" config is appropriate <br>

### make a **config file** for my case and run it: 
```bash
cd "$ALIGNDIR"

# write config file
cat > locarna_scnv.cfg << 'EOF'
max-diff-am: 25
max-diff:    60
min-prob:    0.01
plfold-span: 100
indel:       -50
indel-open:  -750
threads:     2
alifold-consensus-dp
stockholm
EOF
```
explanation: 

```bash
# run
mlocarna \
  --configure locarna_scnv.cfg \
  --tgtdir scnv_min4_locarna \
  "$MUTSEQ/scnv_mutseqs.min4.fasta"
```
started at 12:43 pm <br>
run took a while and srun ended at 3:26.. <br>
restart, but skip making the dot brackets because that step was already done in inputs dir <br>
```bash
srun --pty -t 08:00:00 -c 4 --mem=32G \
  --chdir=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT /bin/bash -l
```

```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
SCRIPTS="$ROOT/scripts"
MUTSEQ="$ROOT/work/06_mutseqs/260403"
ALIGNDIR="$MUTSEQ/07_align"
mkdir -p "$ALIGNDIR"
cd "$ALIGNDIR"

conda activate viromelib


mlocarna \
  --configure locarna_scnv.cfg \
  --skip-pp \
  --tgtdir scnv_min4_locarna \
  "$MUTSEQ/scnv_mutseqs.min4.fasta"

```
`--skip-pp` – won't redo the dot plots 
started 3:29 <br>

### took way too long! 
### next option, filter by uniqueness and bring alignment to ~200 unique seqs 
use `cd-hit`
```bash
conda install -c bioconda cd-hit -y
```

```bash
cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs.min4.fasta" \
  -o "$MUTSEQ/scnv_mutseqs.min4.cdhit95.fasta" \
  -c 0.95 \
  -n 8 \
  -T 4 \
  -M 16000 \
  -d 0

grep -c "^>" "$MUTSEQ/scnv_mutseqs.min4.cdhit95.fasta"
```
```bash
# try 94% identity
cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs.min4.fasta" \
  -o "$MUTSEQ/scnv_mutseqs.min4.cdhit94.fasta" \
  -c 0.94 \
  -n 8 \
  -T 4 \
  -M 16000 \
  -d 0

grep -c "^>" "$MUTSEQ/scnv_mutseqs.min4.cdhit94.fasta"
```

Now, with 245 clusters, use this as the seed for mLocaRNA
```bash
mlocarna \
  --configure locarna_scnv.cfg \
  --tgtdir scnv_min4_locarna_245 \
  "$MUTSEQ/scnv_mutseqs.min4.cdhit94.fasta"
```

change config : 
```bash
cat > locarna_scnv.cfg << 'EOF'
max-diff-am: 25
max-diff:    60
min-prob:    0.01
indel:       -50
indel-opening: -750
threads:     2
alifold-consensus-dp
stockholm
EOF
```
```bash
# check fasta formatting 
head "$MUTSEQ/scnv_mutseqs.min4.cdhit94.fasta"
>seq000001_nmuts5_38,63,80,81,85
CCGTCAGCCTGACGTTAACTGCCACTTTGGTTGGAGTGTGATCAACCGTGTGTGGTGCGAGCATATCGGCCC
>seq000002_nmuts5_38,53,63,66,77
CCGTCAGCGTGACGTTAACTGCCCCTTTGGTTGTAGGGTGATCAACCTTGCCTGGGGCGAGCATATCGGCCC
>seq000003_nmuts5_56,57,79,83,86
CCGTCAGCTTGACGTTAACTGCCACTCGGGTTGAAGTGTGATCAACCGTCCCTTGGACGAGCATATCGGCCC
>seq000004_nmuts5_63,80,81,85,86
CCGTCAGCTTGACGTTAACTGCCACTTTGGTTGTAGTGTGATCAACCGTGAGTGGACCGAGCATATCGGCCC
>seq000005_nmuts5_35,57,75,82,88
CCGTCCGCTTGACGTTAACTGCCACTTGGGTTGAAGTGTGATCAATCGTGCCAGGGGCCAGCATATCGGCCC
```
modify fasta file: 
```bash
sed 's/>\(seq[0-9]*_nmuts[0-9]*\)_.*/>\1/' \
  "$MUTSEQ/scnv_mutseqs.min4.cdhit94.fasta" \
  > "$MUTSEQ/scnv_mutseqs.min4.cdhit94.clean.fasta"

# verify
grep "^>" "$MUTSEQ/scnv_mutseqs.min4.workflows/notescdhit94.clean.fasta" | head -5
```

```bash
mlocarna \
  --configure locarna_scnv.cfg \
  --verbose \
  --tgtdir scnv_min4_locarna_245 \
  "$MUTSEQ/scnv_mutseqs.min4.cdhit94.clean.fasta"
  ```

output 
```
alifold            .(((((...)))))......(.(((...((((((.......))))))))).).(((.(((.....))).))) (-24.46 = -15.14 +  -9.31)
```
install Rscape 
```bash
conda install -c bioconda rscape -y
```

try to get more seq. diversity with another cd hit clustering: 
```bash
cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs.min4.fasta" \
  -o "$MUTSEQ/scnv_mutseqs.min4.cdhit97.fasta" \
  -c 0.97 \
  -n 8 \
  -T 4 \
  -M 16000 \
  -d 0

grep -c "^>" "$MUTSEQ/scnv_mutseqs.min4.cdhit97.fasta"
```
```bash
sed 's/>\(seq[0-9]*_nmuts[0-9]*\)_.*/>\1/' \
  "$MUTSEQ/scnv_mutseqs.min4.cdhit97.fasta" \
  > "$MUTSEQ/scnv_mutseqs.min4.cdhit97.clean.fasta"

```

```bash
mlocarna \
  --configure locarna_scnv.cfg \
  --verbose \
  --tgtdir scnv_min4_locarna_245 \
  "$MUTSEQ/scnv_mutseqs.min4.cdhit97.clean.fasta"
```

## R-scape
```bash
cd "$ALIGNDIR"
mkdir -p rscape_out

R-scape \
  --outdir rscape_out \
  scnv_min4_locarna_245/results/result.stk

ls rscape_out/
```

```bash 
cat rscape_out/result_1.cov
```
output
```
(viromelib) 08:30:13 jgg2144@c0703:07_align# cat rscape_out/result_1.cov
#
# Method Target_E-val [cov_min,cov_max] [FP | TP True Found | Sen PPV F] 
# GTp    0.05           [-3.86,28.98]    [0 | 0 20 0 | 0.00 0.00 0.00] 
#
#-------------------------------------------------------------------------------------------------------
no significant pairs
(viromelib) 08:30:55 jgg2144@c0703:07_align# 
```
look at ranking of bps for closest to significance 
```bash
cat rscape_out/result_1.power
```

# Refine mutants: 2 muts 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
SCRIPTS="$ROOT/scripts"
TRIM="$ROOT/work/02_trim/260327" # for date 
ALIGN="$ROOT/work/04_align/260327" # for date 
REFDIR="$ROOT/data/refseqs"

MUTSEQ="$ROOT/work/06_mutseqs/260403"
ALIGNDIR="$MUTSEQ/07_align"
mkdir -p "$ALIGNDIR"
cd "$ALIGNDIR"

python "$SCRIPTS/extract_mutseqs.py" \
  --bam-x1  "$ALIGN/260325-r0415-S-X-1-A3408_R1.NOLEAD.sorted.bam" \
  --bam-x2  "$ALIGN/260325-r0415-S-X-2-A3477_R1.NOLEAD.sorted.bam" \
  --ref     "$REFDIR/scnv_xr2_wt.fa" \
  --pos-min 30 \
  --pos-max 101 \
  --out     "$MUTSEQ/scnv_mutseqs" \
  --min-muts 2 \ 
  --max-seqs 10000

```
`min-muts` just sets a floor , 
To get only exactly 2-mutation sequences, filter after:
```bash
# keep only nmuts2 sequences from the output
grep -A1 "nmuts2" "$MUTSEQ/scnv_mutseqs.min2.fasta" | \
  grep -v "^--$" > "$MUTSEQ/scnv_mutseqs.exact2.fasta"

grep -c "^>" "$MUTSEQ/scnv_mutseqs.exact2.fasta"
```
## now do cd hit 
```bash
cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs.exact2.fasta" \
  -o "$MUTSEQ/scnv_mutseqs.exact2.cdhit97.fasta" \
  -c 0.97 \
  -n 8 -T 4 -M 16000 -d 0

grep -c "^>" "$MUTSEQ/scnv_mutseqs.exact2.cdhit97.fasta"
```

## try mLocaRNA without the cluster reduction
```bash
# clean headers 
sed 's/>\(seq[0-9]*_nmuts[0-9]*\)_.*/>\1/' \
  "$MUTSEQ/scnv_mutseqs.exact2.fasta" \
  > "$MUTSEQ/scnv_mutseqs.exact2.clean.fasta"

# run
mlocarna \
  --configure "$ALIGNDIR/locarna_scnv.cfg" \
  --verbose \
  --tgtdir scnv_exact2_locarna \
  "$MUTSEQ/scnv_mutseqs.exact2.clean.fasta"
```

too many seqs. took way too long. do the 98 seqs 

```bash
# Ctrl+C to kill mlocarna

# clean headers on the 98-seq file
sed 's/>\(seq[0-9]*_nmuts[0-9]*\)_.*/>\1/' \
  "$MUTSEQ/scnv_mutseqs.exact2.cdhit97.fasta" \
  > "$MUTSEQ/scnv_mutseqs.exact2.cdhit97.clean.fasta"

mlocarna \
  --configure "$ALIGNDIR/locarna_scnv.cfg" \
  --verbose \
  --tgtdir scnv_exact2_98seqs \
  "$MUTSEQ/scnv_mutseqs.exact2.cdhit97.clean.fasta"
```
### `R-scape`
```bash
# once mlocarna finishes, check the consensus structure
grep "SS_cons" scnv_exact2_98seqs/results/result.stk

# run R-scape within $ALIGNDIR
mkdir -p rscape_exact2
R-scape \
  --outdir rscape_exact2 \
  scnv_exact2_98seqs/results/result.stk

# check results
cat rscape_exact2/*.cov
cat rscape_exact2/*.power
```

# Mut burden plot 
### Run the `bam_wt_mut_burden.py` script on beny and scnv samples 
### Apr 7, 2026 
```bash
BURDEN="$ROOT/work/08_wt_freq/260327"
mkdir -p "$BURDEN"

for BAM in "$ALIGN"/*.sorted.bam; do
  BASE=$(basename "$BAM" .sorted.bam)
  [[ "$BASE" == *-R-* && "$BASE" != *.HASLEAD ]] && continue
  [[ "$BASE" == *-X-* && "$BASE" != *.NOLEAD  ]] && continue
  [[ "$BASE" == *-S-* ]] && REF="$REFDIR/scnv_xr2_wt.fa" && PMIN=30 && PMAX=101
  [[ "$BASE" == *-B-* ]] && REF="$REFDIR/beny_xr_wt.fa"  && PMIN=30 && PMAX=101

  echo "Processing $BASE"
  python "$SCRIPTS/bam_wt_mut_burden.py" \
    "$BAM" "$REF" "$PMIN" "$PMAX" \
    "$BURDEN/$BASE"
done
```
```
bam_wt_mut_burden:  per READ  — how many mutations does this read carry overall?
wtnorm_fc:          per POSITION — how does mutation frequency at this spot change?
```

### Make metadata table to run mut burden in R 
```bash
# create sample table
cat > "$BURDEN/sample_table_edit.csv" << 'EOF'
sample,virus,condition,replicate
260325-r0415-S-X-1-A3408_R1.NOLEAD,SCNV,X,1
260325-r0415-S-X-2-A3477_R1.NOLEAD,SCNV,X,2
260325-r0415-S-R-1-A3449_R1.HASLEAD,SCNV,R,1
260325-r0415-S-R-2-A3405_R1.HASLEAD,SCNV,R,2
260325-r0414-B-X-1-A3471_R1.NOLEAD,BENY,X,1
260325-r0414-B-X-2-A3454_R1.NOLEAD,BENY,X,2
260325-r0414-B-R-1-A3474_R1.HASLEAD,BENY,R,1
260325-r0414-B-R-2-A3470_R1.HASLEAD,BENY,R,2
EOF
```
```bash

# scp to local
# downloaded to ~/Desktop/viromelib/htpmut/work/08_wt_freq/260327

# run R locally

```

# `mLocaRNA` - SCNV core 
### Apr 7, 2026

shorter sequence = faster mLocaRNA ... potentially stronger covariation if focusing on the structural core <br>
Re-extract mut counts with the tighter window: <br>
```bash 
python "$SCRIPTS/extract_mutseqs.py" \
  --bam-x1  "$ALIGN/260325-r0415-S-X-1-A3408_R1.NOLEAD.sorted.bam" \
  --bam-x2  "$ALIGN/260325-r0415-S-X-2-A3477_R1.NOLEAD.sorted.bam" \
  --ref     "$REFDIR/scnv_xr2_wt.fa" \
  --pos-min 30 \
  --pos-max 80 \
  --out     "$MUTSEQ/scnv_mutseqs_core" \
  --min-muts 2 \
  --max-seqs 10000
```
high amount - `skipped non-spanning`:1766704 <br>
new script spanning check is replaced with an overlap check
```
# OLD: strict full span required
if read.reference_start > (pos_min - 1) or read.reference_end < pos_max:
    skipped_span += 1
    continue

# NEW: just needs to overlap the window at all
if read.reference_end <= (pos_min - 1) or read.reference_start >= pos_max:
    skipped_no_overlap += 1
    continue
```
Uncovered positions within the window stay as reference bases — so two reads with the same mutations but different coverage extents will produce the same sequence string and deduplicate correctly.

```bash
python "$SCRIPTS/extract_mutseqs_partial.py" \
  --bam-x1  "$ALIGN/260325-r0415-S-X-1-A3408_R1.NOLEAD.sorted.bam" \
  --bam-x2  "$ALIGN/260325-r0415-S-X-2-A3477_R1.NOLEAD.sorted.bam" \
  --ref     "$REFDIR/scnv_xr2_wt.fa" \
  --pos-min 30 \
  --pos-max 80 \
  --out     "$MUTSEQ/scnv_mutseqs_core_partial" \
  --min-muts 2 3 \
  --max-seqs 10000
```

```bash
# try different thresholds on the partial output
cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs_core_partial.min2.fasta" \
  -o "$MUTSEQ/scnv_mutseqs_core_partial.cdhit94.fasta" \
  -c 0.94 -n 8 -T 4 -M 16000 -d 0
grep -c "^>" "$MUTSEQ/scnv_mutseqs_core_partial.cdhit94.fasta"

cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs_core_partial.min2.fasta" \
  -o "$MUTSEQ/scnv_mutseqs_core_partial.cdhit96.fasta" \
  -c 0.96 -n 8 -T 4 -M 16000 -d 0
grep -c "^>" "$MUTSEQ/scnv_mutseqs_core_partial.cdhit96.fasta"

cd-hit-est \
  -i "$MUTSEQ/scnv_mutseqs_core_partial.min2.fasta" \
  -o "$MUTSEQ/scnv_mutseqs_core_partial.cdhit98.fasta" \
  -c 0.98 -n 8 -T 4 -M 16000 -d 0
grep -c "^>" "$MUTSEQ/scnv_mutseqs_core_partial.cdhit98.fasta"
```

output to move forward with: 
```
Command: cd-hit-est -i
         /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs_core_partial.min2.fasta
         -o
         /groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_mutseqs/260403/scnv_mutseqs_core_partial.cdhit94.fasta
         -c 0.94 -n 8 -T 4 -M 16000 -d 0

Started: Tue Apr  7 17:12:10 2026
================================================================
                            Output                              
----------------------------------------------------------------
total seq: 10000
longest and shortest : 51 and 51
Total letters: 510000
Sequences have been sorted

Approximated minimal memory consumption:
Sequence        : 1M
Buffer          : 4 X 17M = 70M
Table           : 2 X 1M = 2M
Miscellaneous   : 0M
Total           : 75M

Table limit with the given memory limit:
Max number of representatives: 4000000
Max number of word counting entries: 1990567023

# comparing sequences from          0  to       1666
.---------- new table with      287 representatives
# comparing sequences from       1666  to       3055
99.9%---------- new table with       42 representatives
# comparing sequences from       3055  to       4212
95.0%---------- new table with       26 representatives
# comparing sequences from       4212  to       5176
95.0%---------- new table with       25 representatives
# comparing sequences from       5176  to       5980
72.0%---------- new table with        9 representatives
# comparing sequences from       5980  to       6650
---------- new table with       10 representatives
# comparing sequences from       6650  to       7208
---------- new table with        7 representatives
# comparing sequences from       7208  to       7673
---------- new table with        1 representatives
# comparing sequences from       7673  to       8060
---------- new table with        0 representatives
# comparing sequences from       8060  to       8383
---------- new table with        0 representatives
# comparing sequences from       8383  to       8652
---------- new table with        0 representatives
# comparing sequences from       8652  to       8876
---------- new table with        0 representatives
# comparing sequences from       8876  to       9063
---------- new table with        0 representatives
# comparing sequences from       9063  to      10000
---------- new table with        0 representatives

    10000  finished        407  clusters

```

### clean up 94 file: 
 ```bash
sed 's/>\(seq[0-9]*_nmuts[0-9]*\)_.*/>\1/' \
  "$MUTSEQ/scnv_mutseqs_core_partial.cdhit94.fasta" \
  > "$MUTSEQ/scnv_mutseqs_core_partial.cdhit94.clean.fasta"

mlocarna \
  --configure "$ALIGNDIR/locarna_scnv.cfg" \
  --verbose \
  --tgtdir "$ALIGNDIR/scnv_core_partial_locarna" \
  "$MUTSEQ/scnv_mutseqs_core_partial.cdhit94.clean.fasta"

 ```

 Note - make sure position directions in script are saying what I want ("look for mutations in this position range"), not "only keep reads that start and end at these sites" 

### `R-scape` with ≥2 mut , `CD-HIT`-processed, `LocARNA` aligned reads. 

```bash
mkdir -p "$ALIGNDIR/rscape_core_partial"

R-scape \
  --outdir "$ALIGNDIR/rscape_core_partial" \
  "$ALIGNDIR/scnv_core_partial_locarna/results/result.stk"

cat "$ALIGNDIR/rscape_core_partial/"*.cov
cat "$ALIGNDIR/rscape_core_partial/"*.power
```
power analysis is now much better than before — pairs 23-50, 29-47, 31-45, 32-44 all have power >0.85 = enough sequences and substitutions to detect covariation ... But observed 0 covarying pairs out of expected 10.6. <br>


With 407 sequences, 51nt window, and power >0.85 at the key pairs. <br>
The absence of covariation is real signal, not a power problem.

likely due to... complex structure. The scnv xrRNA likely has complex interactions rather than simple nested stems <br>
not just pairing


