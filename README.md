# HTP MUT 
## analysis of depleted mutant xrRNA sequences following exonuclease treatment 

## set up project dir
**data stored in** `/groups/as6282_gp/data_bkup/jgg2144/20260121_Aviti/`
**project directory stored in** `/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT`

`ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/` <br>

`"$ROOT"/data` contains: <br>
- `"$ROOT"/data/refseqs`

`"$ROOT"/work` contains: <br>
- `"$ROOT"/work/01_raw`: symlinked raw fastq files
- `"$ROOT"/work/02_trim`: trimmed reads using cutadapt
- `"$ROOT"/work/03_qc`: post-trim qc using fastqc multiqc 
- `"$ROOT"/work/04_align`: alignment to the wild type using bowtie2

`"$ROOT"/scripts` contains: <br>


## Interactive job (`srun --pty`)
```bash
srun --pty -t 06:00:00 -c 4 --mem=16G \
--chdir=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT \
/bin/bash -l
```

# 01_raw
```bash
export ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
export RAW="$ROOT/work/01_raw"
mkdir -p "$RAW"
cd "$RAW"
```
### symlink 
```bash
export SRC=/groups/as6282_gp/data_bkup/jgg2144/20260121_Aviti
ls $SRC/*.fastq.gz
```

### had previously set up a conda env, for `viromelib` project, will activate: 
**(done before):** <br>
```bash
conda create -y -n viromelib -c conda-forge -c bioconda \
  python=3.11 cutadapt fastp bowtie2 samtools fastqc multiqc seqkit
```
```bash
# in viromelib conda env 
conda install -c conda-forge ncbi-datasets-cli -y
conda install -c bioconda bedtools -y
```
**now:** <br>
```bash
conda activate viromelib
```

# 02_trim
```bash 
TRIM=$ROOT/work/02_trim
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
cd "$TRIM"

for SAMPLE in "$RAW"/*_R1.fastq.gz; do
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

outputs: <br>
```
Done           00:00:09     5,323,951 reads @   1.8 µs/read;  34.06 M reads/minute
Done           00:00:10     5,824,610 reads @   1.8 µs/read;  33.91 M reads/minute
Done           00:00:08     4,850,366 reads @   1.8 µs/read;  32.45 M reads/minute
Done           00:00:10     6,626,167 reads @   1.5 µs/read;  39.25 M reads/minute
Done           00:00:09     6,031,584 reads @   1.5 µs/read;  39.37 M reads/minute
Done           00:00:10     6,940,111 reads @   1.5 µs/read;  40.16 M reads/minute
Done           00:00:05     2,834,083 reads @   1.9 µs/read;  32.36 M reads/minute
Done           00:00:07     3,833,112 reads @   1.8 µs/read;  32.56 M reads/minute
Done           00:00:08     4,496,671 reads @   1.9 µs/read;  31.94 M reads/minute
Done           00:00:09     5,827,147 reads @   1.7 µs/read;  35.05 M reads/minute
Done           00:00:09     5,542,556 reads @   1.7 µs/read;  34.29 M reads/minute
Done           00:00:10     6,368,859 reads @   1.7 µs/read;  35.02 M reads/minute
```

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
nano "$SCRIPTS/leader_summary.sh"
bash "$SCRIPTS/leader_summary.sh"
```
Loops over every *.HASLEAD.fastq.gz file in `02_trim` ; Finds the matching .NOLEAD.fastq.gz for each sample; Counts reads (lines / 4); Computes total reads per sample; Computes fraction of reads retaining the leader (HASLEAD / Total); Writes all of this to leader_summary.tsv

# 03_qc
```bash
TRIM=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/02_trim
QC=$ROOT/work/03_qc
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
```bash
# scp to local and view report
# on local machine; 
scp -r jgg2144@hpc.c2b2.columbia.edu:/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/03_qc/multiqc_report.html ~/Downloads/
```
# 04_align
```bash 
mkdir -p "$ROOT/data/refseqs"
cd "$ROOT/data/refseqs"
nano refseqs/zikv_xr2_wt.fa
nano refseqs/denv_xr2_wt.fa
```
for now, both references include universal 5' end, but not the universal 3' end. <br>
**`zikv_xr2_wt.fa`** <br>
`GGATTAATATAATCTGGGAAACCAAGCTCATAGTCAGGCCGAGAACGCCATGGCACGGAAGAAGCCATGCTGCCTGTGAGCCCCTCAGAGGACACTGAGTCAAAAAACCCCAC` <br>
**`denv_xr2_wt.fa`**<br>
`GGATTAATATAATCCAAGGACGTTAAAAGAAGTCAGGCCATCATAAATGCCATAGCTGGAGTAAACTATGCAGCCTGTAGCTCCACCTGAGAAGGTGTAAAAAATCCGGGAGG`

– notes - <br>
-will use `--very-sensitive-local` so that Bowtie2 can align part of the read, even if it isn't a perfect match (contains muts.)<br>
-recall unmatched ends will be soft-clipped <br>
–recall that ligation oligo added max 8 bases, which will be soft-clipped in either case <br>

**index** - create the binary lookup tables Bowtie2 uses for fast matching: <br>
```bash
cd "$ROOT/data/refseqs"
bowtie2-build zikv_xr2_wt.fa zikv_xr2_wt
bowtie2-build denv_xr2_wt.fa denv_xr2_wt
```
now, navigate to align dir: 
```bash
ALIGN="$ROOT/work/04_align"
TRIM="$ROOT/work/02_trim"
mkdir -p "$ALIGN"
cd "$ALIGN"
```
### try align with one sample pair first: 
```bash
SAMPLE=260120-r0399-ZX-01-A3466_R1
REF="$ROOT/data/refseqs/zikv_xr2_wt"
```

```bash 
bowtie2 -p 8 \
  --very-sensitive-local \
  -x "$REF" \
  -U "$ROOT/work/02_trim/${SAMPLE}.insert.fastq.gz" \
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
example output: <br>
```99.09% overall alignment rate
...skipping...
3673782 reads; of these:
  3673782 (100.00%) were unpaired; of these:
    33262 (0.91%) aligned 0 times
    3638607 (99.04%) aligned exactly 1 time
    1913 (0.05%) aligned >1 times
99.09% overall alignment rate
```
do for all 12 samples independently. <br>

### or, can loop alignment, sorting, indexing, flagstat: <br>
```bash
TRIM="$ROOT/work/02_trim"
REFDIR="$ROOT/data/refseqs"

for READS in "$TRIM"/*.insert.fastq.gz; do
  SAMPLE=$(basename "$READS" .insert.fastq.gz)

  # choose reference based on virus
  if [[ "$SAMPLE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt"
  elif [[ "$SAMPLE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt"
  else
    echo "Unknown virus for $SAMPLE — skipping"
    continue
  fi

  echo "Aligning $SAMPLE"

  bowtie2 -p 8 \
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

check overall stat:
```bash
grep "overall alignment rate" *.bowtie2.log
```
example output: 
```
260120-r0399-ZR-01-A3483_R1.bowtie2.log:99.86% overall alignment rate
260120-r0399-ZR-02-A3456_R1.bowtie2.log:99.87% overall alignment rate
260120-r0399-ZR-03-A3451_R1.bowtie2.log:99.86% overall alignment rate
260120-r0399-ZX-01-A3466_R1.bowtie2.log:99.09% overall alignment rate
260120-r0399-ZX-02-A3452_R1.bowtie2.log:99.18% overall alignment rate
260120-r0399-ZX-03-A3450_R1.bowtie2.log:99.01% overall alignment rate
260120-r0400-DR-01-A3413_R1.bowtie2.log:99.92% overall alignment rate
260120-r0400-DR-02-A3419_R1.bowtie2.log:99.91% overall alignment rate
260120-r0400-DR-03-A3485_R1.bowtie2.log:99.89% overall alignment rate
260120-r0400-DX-01-A3429_R1.bowtie2.log:99.74% overall alignment rate
260120-r0400-DX-02-A3499_R1.bowtie2.log:99.71% overall alignment rate
260120-r0400-DX-03-A3484_R1.bowtie2.log:99.66% overall alignment rate
```
# 04_align for leader-split reads 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT

TRIM="$ROOT/work/02_trim"        # trimmed FASTQs (HASLEAD + NOLEAD)
REFDIR="$ROOT/data/refseqs"      # bowtie2 indexes
ALIGN="$ROOT/work/04_align"

mkdir -p "$ALIGN"
cd "$ALIGN"
```
alignment loop for both leader-containing and no-lead samples after the leader-split step: <br>
``` bash
for READS in "$TRIM"/*.HASLEAD.fastq.gz "$TRIM"/*.NOLEAD.fastq.gz; do
  [[ -e "$READS" ]] || continue

  SAMPLE=$(basename "$READS" .fastq.gz)

  # Choose reference by virus label in filename
  if [[ "$SAMPLE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt"
  elif [[ "$SAMPLE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt"
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
Note: current settings include NOLEAD for R and HASLEAD for X, which is not necessary to do and should have very minimal reads <br>
```
260120-r0399-ZR-01-A3483_R1.HASLEAD.bowtie2.log:99.92% overall alignment rate
260120-r0399-ZR-01-A3483_R1.NOLEAD.bowtie2.log:98.18% overall alignment rate
260120-r0399-ZR-01-A3483_R1.bowtie2.log:99.86% overall alignment rate
260120-r0399-ZR-02-A3456_R1.HASLEAD.bowtie2.log:99.92% overall alignment rate
260120-r0399-ZR-02-A3456_R1.NOLEAD.bowtie2.log:98.09% overall alignment rate
260120-r0399-ZR-02-A3456_R1.bowtie2.log:99.87% overall alignment rate
260120-r0399-ZR-03-A3451_R1.HASLEAD.bowtie2.log:99.91% overall alignment rate
260120-r0399-ZR-03-A3451_R1.NOLEAD.bowtie2.log:98.22% overall alignment rate
260120-r0399-ZR-03-A3451_R1.bowtie2.log:99.86% overall alignment rate
260120-r0399-ZX-01-A3466_R1.HASLEAD.bowtie2.log:99.07% overall alignment rate
260120-r0399-ZX-01-A3466_R1.NOLEAD.bowtie2.log:99.10% overall alignment rate
260120-r0399-ZX-01-A3466_R1.bowtie2.log:99.09% overall alignment rate
260120-r0399-ZX-02-A3452_R1.HASLEAD.bowtie2.log:98.98% overall alignment rate
260120-r0399-ZX-02-A3452_R1.NOLEAD.bowtie2.log:99.18% overall alignment rate
260120-r0399-ZX-02-A3452_R1.bowtie2.log:99.18% overall alignment rate
260120-r0399-ZX-03-A3450_R1.HASLEAD.bowtie2.log:98.45% overall alignment rate
260120-r0399-ZX-03-A3450_R1.NOLEAD.bowtie2.log:99.03% overall alignment rate
260120-r0399-ZX-03-A3450_R1.bowtie2.log:99.01% overall alignment rate
260120-r0400-DR-01-A3413_R1.HASLEAD.bowtie2.log:99.95% overall alignment rate
260120-r0400-DR-01-A3413_R1.NOLEAD.bowtie2.log:98.54% overall alignment rate
260120-r0400-DR-01-A3413_R1.bowtie2.log:99.92% overall alignment rate
260120-r0400-DR-02-A3419_R1.HASLEAD.bowtie2.log:99.95% overall alignment rate
260120-r0400-DR-02-A3419_R1.NOLEAD.bowtie2.log:98.42% overall alignment rate
260120-r0400-DR-02-A3419_R1.bowtie2.log:99.91% overall alignment rate
260120-r0400-DR-03-A3485_R1.HASLEAD.bowtie2.log:99.96% overall alignment rate
260120-r0400-DR-03-A3485_R1.NOLEAD.bowtie2.log:98.25% overall alignment rate
260120-r0400-DR-03-A3485_R1.bowtie2.log:99.89% overall alignment rate
260120-r0400-DX-01-A3429_R1.HASLEAD.bowtie2.log:99.73% overall alignment rate
260120-r0400-DX-01-A3429_R1.NOLEAD.bowtie2.log:99.74% overall alignment rate
260120-r0400-DX-01-A3429_R1.bowtie2.log:99.74% overall alignment rate
260120-r0400-DX-02-A3499_R1.HASLEAD.bowtie2.log:99.75% overall alignment rate
260120-r0400-DX-02-A3499_R1.NOLEAD.bowtie2.log:99.71% overall alignment rate
260120-r0400-DX-02-A3499_R1.bowtie2.log:99.71% overall alignment rate
260120-r0400-DX-03-A3484_R1.HASLEAD.bowtie2.log:99.73% overall alignment rate
260120-r0400-DX-03-A3484_R1.NOLEAD.bowtie2.log:99.66% overall alignment rate
260120-r0400-DX-03-A3484_R1.bowtie2.log:99.66% overall alignment rate
```
# 05_readcount (with updates) 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT

ALIGN="$ROOT/work/04_align"          # where *.sorted.bam live
OUT="$ROOT/work/05_readcount"        # output directory
REFDIR="$ROOT/data/refseqs"          # fasta references (not bowtie2 indexes)

mkdir -p "$OUT"
```
### loop: <br>
```bash
for BAM in "$ALIGN"/*.sorted.bam; do
  [[ -e "$BAM" ]] || continue
  BASE=$(basename "$BAM" .sorted.bam)

  # Choose FASTA reference based on virus label in filename
  if [[ "$BASE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt.fa"
  elif [[ "$BASE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt.fa"
  else
    echo "Unknown sample: $BASE — skipping"
    continue
  fi

  # Determine contig name + length from BAM index stats
  IDXLINE=$(samtools idxstats "$BAM" | head -n 1)
  CONTIG=$(echo "$IDXLINE" | cut -f1)
  LEN=$(echo "$IDXLINE" | cut -f2)
  REGION="${CONTIG}:1-${LEN}"

  echo "bam-readcount: $BASE on $REGION"

  bam-readcount \
    -f "$REF" \
    -d 10000000 \
    -q 20 \
    -b 20 \
    "$BAM" "$REGION" \
    > "$OUT/${BASE}.readcount.txt"
done
```
Notes: <br>
- takes about 4 min per aligned sample of around 4 mill reads. <br>
- should probably change naming convention to be more parsable and not rely on just one letter Z or D for example.. <br>
- `-f "$REF"` is the FASTA reference sequence 
- `-d 10000000` - Maximum depth per position (set higher than my highest read depth); Default is often too low for ultra-deep libraries.
- `-q 20` - Minimum mapping quality (MAPQ); Filters reads that align ambiguously. Ref is small size, and mapping rate is high, so MAPQ should be pretty high
- `-b 20` - Minimum base quality - prevents low-quality sequencing errors from inflating mutation counts.
- `"$BAM" "$REGION"` - Run only over the single contig region (1 to full length) Faster than scanning whole BAM blindly. Bc of idxstats, region is auto-detected from the BAM 

### same loop but only for HASLEAD for No-Xrn1, NOLEAD for Xrn1: 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
OUT="$ROOT/work/05_readcount"
REFDIR="$ROOT/data/refseqs"

mkdir -p "$OUT"

for BAM in "$ALIGN"/*.sorted.bam; do
  [[ -e "$BAM" ]] || continue
  BASE=$(basename "$BAM" .sorted.bam)

  # Keep only:
  #   R samples -> HASLEAD
  #   X samples -> NOLEAD
  if [[ "$BASE" == *-R-* ]]; then
    [[ "$BASE" == *.HASLEAD ]] || continue
  elif [[ "$BASE" == *-X-* ]]; then
    [[ "$BASE" == *.NOLEAD ]] || continue
  else
    echo "Skipping (not R or X): $BASE"
    continue
  fi

  # Reference fasta based on virus
  if [[ "$BASE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt.fa"
  elif [[ "$BASE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt.fa"
  else
    echo "Unknown virus for $BASE — skipping"
    continue
  fi

  # Determine contig name + length from BAM index stats
  IDXLINE=$(samtools idxstats "$BAM" | head -n 1)
  CONTIG=$(echo "$IDXLINE" | cut -f1)
  LEN=$(echo "$IDXLINE" | cut -f2)
  REGION="${CONTIG}:1-${LEN}"

  echo "bam-readcount: $BASE on $REGION"

  bam-readcount \
    -f "$REF" \
    -d 10000000 \
    -q 20 \
    -b 20 \
    "$BAM" "$REGION" \
    > "$OUT/${BASE}.readcount.txt"
done
```
`* -R-*` samples: only runs if filename ends with .HASLEAD; 
`* -X-*` samples: only runs if filename ends with .NOLEAD

# 06_counts

## transform BAM readcounts into mutation features for frequency analysis
```
bam-readcount outputs
        ↓
mutation features (pos_ref>alt)
        ↓
counts matrix
        ↓
DESeq2 statistics
        ↓
plot
```
```
06_counts/
 ── alt_counts.tsv     # matrix for DESeq2
 ── sample_table.tsv   # metadata
 ```
 
create output dir: 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
mkdir -p $ROOT/work/06_counts
```
run python script **readcount_to_matrix.py**: <br>
```bash
export READCOUNT_DIR=$ROOT/work/05_readcount
export OUTDIR=$ROOT/work/06_counts

python $ROOT/scripts/readcount_to_matrix.py
```

output files: <br>
`alt_counts.tsv`: Rows = mutation features
example: 
```
feature        ZR1   ZR2   ZX1 ...
31_A>G         523   510   22
31_A>C         490   501   430
...
```
Each row = one substituion mutation at given position <br>

`sample_table.tsv`
```
sample     virus   condition   replicate
ZR-01      ZIKV    R           01
ZX-01      ZIKV    X           01
...
```
tells DESeq2 what to compare 

# 07_DESeq2
scp to local and open R. 
```bash
scp -r jgg2144@hpc.c2b2.columbia.edu:/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/06_counts/"*.tsv" work/06_counts
```
manually edited the data files so that only "HASLEAD" for No-Xrn1 and "NOLEAD" for Xrn1-treated were included for each replicate. 

find DESeq2 and plotting information in `/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/workflows/DESeq2/mut_count_plot_2.rmd`

# anchor to invariant WT region , and have DESeq do ratio of ratios for each condition relative to WT. 
### relies on a scripts named `bam_wt_mut_burden`, which: 
- restricts to the mutated/degenerate base window (for ZIKV, 29-107, for DENV, 28-104)
- if a read base is not equal to the reference base, call as mutation 
- qc filters: base quality ≥ 20, mapping quality ≥ 20
- count mutations per read 
- get WT frequency, `wt_reads`= number of reads with 0 mutations; `wt_fraction`= `wt_reads`/`total_reads` - giving the WT haplotype abundance
- also get the mutation burden: frequency of n=1,2,3.. mutations (get an idea of how many mutations coexist on same molecule)
- results in a true WT reference, direct measure of exact WT on the SEQUENCE level, not BASE level.
- can now normalize the per BASE mutation frequency to the WT frequency.
- log2FC will now be:
```
log2FC = (mutant/WT)_X vs (mutant/WT)_R
```
- sets WT as baseline, 0, so now when plotting, each dot means: how much this mutation changed in X vs R relative to WT molecules
- can apply the mutation burden data to filter only mutation burden = 1 for single mutants. (will have fewer reads as some signal is lost)




 ```bash
 ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
OUT="$ROOT/work/08_wt_freq"
REFDIR="$ROOT/data/refseqs"
SAMPLETAB="$ROOT/work/06_counts/sample_table_edit.csv" # this csv was manually edited locally, so only the Xrn1-treated samples without leader and the No-Xrn1 samples with the leader were considered for analysis. 

mkdir -p "$OUT"

ZIKV_MIN=29; ZIKV_MAX=107
DENV_MIN=28; DENV_MAX=104

tail -n +2 "$SAMPLETAB" | while IFS=, read -r sample virus condition replicate; do

  BAM="$ALIGN/${sample}.sorted.bam"

  if [[ ! -f "$BAM" ]]; then
    echo "Missing BAM: $BAM" >&2
    continue
  fi

  if [[ "$virus" == "ZIKV" ]]; then
    REF="$REFDIR/zikv_xr2_wt.fa"
    POSMIN=$ZIKV_MIN
    POSMAX=$ZIKV_MAX
  else
    REF="$REFDIR/denv_xr2_wt.fa"
    POSMIN=$DENV_MIN
    POSMAX=$DENV_MAX
  fi

  echo "WT+burden: $sample ($virus $condition)"
  python "$ROOT/scripts/bam_wt_mut_burden.py" \
    "$BAM" "$REF" "$POSMIN" "$POSMAX" \
    "$OUT/$sample"

done
```

fix header of `HTPMUT/work/08_wt_freq` - issue in making the table.. <br>
```
OUT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/08_wt_freq

> "$OUT/WT_summary_all.tsv"

first=1
for f in "$OUT"/*.summary.tsv; do
  if [[ $first -eq 1 ]]; then
    cat "$f" > "$OUT/WT_summary_all.tsv"
    first=0
  else
    tail -n +2 "$f" >> "$OUT/WT_summary_all.tsv"
  fi
done

head -n 5 "$OUT/WT_summary_all.tsv"
wc -l "$OUT/WT_summary_all.tsv"
```

```
# steps: 

Plot mutation burden distribution 

Run DESeq2 with WT offset on all reads

Generate single-mutant matrix

Run DESeq2 on single mutants only

Compare bubble plots
```
### **Follow plotting notebook in R notebook** (found in workflows too)

---------- 
# attempts preceding Feb 23, 2026: <br>
# 05_pileup

use samtools to get pileup information (input is the bam file)<br>

use a python script: **pileup_to_counts.py** to clean samtools column into meaninful data <br>

script explained: key function is `mpileup` in samtools <br>
`samtools mpileup -f ref.fa sample.sorted.bam > sample.pileup`<br>
goes through each position in the reference; looks at all reads overlapping that position; reports: reference base, coverage depth, actual bases observed in the reads<br>

column 5 contains symbols for the actual / observed bases: <br>
`.`	match to reference (forward strand)<BR>
`,`	match to reference (reverse strand)<BR>
`A,C,G,T`	mismatch<BR>
`a,c,g,t`	mismatch on reverse strand<BR>
`^`	start of read<BR>
`$`	end of read<BR>
`+2AG`insertion<BR>
`-1T`	deletion <BR>

in order to make sequence logos, clean column 5, so that provides explicit base counts for each base: 
`A_count`
`C_count`
`G_count`
`T_count`

The loop will, at each ref. pos., parse columns, clean the base string, convert matches (., ,) to reference base, count A/C/G/T, store result

### loop through pileup to get counts (depth is limited to 8K): 
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
PILEUP="$ROOT/work/05_pileup"

for BAM in "$ALIGN"/*.sorted.bam; do

  BASENAME=$(basename "$BAM" .sorted.bam)

  # choose reference
  if [[ "$BASENAME" == *-Z* ]]; then
    REF="$ROOT/data/refseqs/zikv_xr2_wt.fa"
  elif [[ "$BASENAME" == *-D* ]]; then
    REF="$ROOT/data/refseqs/denv_xr2_wt.fa"
  else
    echo "Unknown sample: $BASENAME"
    continue
  fi

  echo "Processing $BASENAME"

  # 1. generate pileup
  samtools mpileup -f "$REF" "$BAM" > "$PILEUP/${BASENAME}.pileup"

  # 2. convert pileup to counts
  python "$ROOT/scripts/pileup_to_counts.py" \
    "$PILEUP/${BASENAME}.pileup" \
    "$PILEUP/${BASENAME}.counts.tsv"

done

```
### try pileup without read depth limitation (takes more memory but preserves depth): EDIT: memory hog, doesn't work
```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
PILEUP="$ROOT/work/05_pileup"

mkdir -p "$PILEUP"

for BAM in "$ALIGN"/*.sorted.bam; do
  BASENAME=$(basename "$BAM" .sorted.bam)

  # choose reference
  if [[ "$BASENAME" == *-Z* ]]; then
    REF="$ROOT/data/refseqs/zikv_xr2_wt.fa"
  elif [[ "$BASENAME" == *-D* ]]; then
    REF="$ROOT/data/refseqs/denv_xr2_wt.fa"
  else
    echo "Unknown sample: $BASENAME"
    continue
  fi

  echo "Processing $BASENAME"

  # 1) mpileup: remove depth cap + filter low quality
  samtools mpileup \
    -d 10000000 \
    -q 20 -Q 20 \
    -f "$REF" \
    "$BAM" > "$PILEUP/${BASENAME}.pileup"

  # 2) convert pileup to counts (explicit output filename)
  python "$ROOT/scripts/pileup_to_counts.py" \
    "$PILEUP/${BASENAME}.pileup" \
    "$PILEUP/${BASENAME}.pileup.counts.tsv"
done
```
check that the pileup outputs exist: <br>
```bash
ls $PILEUP | head
```

# 06_logo
### use the per-position base frequencies from `05_pileup` to create input/output sequence logos for each condition
```bash
conda activate viromelib
conda install -c conda-forge logomaker -y # to have for sequence logo making 
```
python script that gathers from the pileup - mutation frequency per position per replicate (with reference WT base column): <br>

```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
PILEUP="$ROOT/work/05_pileup"
SCRIPTS="$ROOT/scripts"

cd "$SCRIPTS"
python mut_freq_with_ref.py
```
-------------
# 05_read count 
## new approach for taking BAM files to get mut. frequencies
```bash
conda install -c bioconda bam-readcount -y
```
produces per-position depth + A/C/G/T counts without the `mpileup` huge “bases string”

```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
ALIGN="$ROOT/work/04_align"
OUT="$ROOT/work/05_readcount"
REFDIR="$ROOT/data/refseqs"

mkdir -p "$OUT"

for BAM in "$ALIGN"/*.sorted.bam; do
  BASE=$(basename "$BAM" .sorted.bam)

  if [[ "$BASE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt.fa"
  elif [[ "$BASE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt.fa"
  else
    echo "Unknown sample: $BASE"
    continue
  fi

  # BAM has one reference contig; grab its name + length
  IDXLINE=$(samtools idxstats "$BAM" | head -n 1)
  CONTIG=$(echo "$IDXLINE" | cut -f1)
  LEN=$(echo "$IDXLINE" | cut -f2)

  REGION="${CONTIG}:1-${LEN}"
 

  echo "bam-readcount: $BASE on $REGION"
  bam-readcount -d 10000000 -f "$REF" -q 20 -b 20 "$BAM" "$REGION" \
    > "$OUT/${BASE}.readcount.txt"
done
```
note - alignment was done with the **.insert.fastq.gz** trimmed files - not the leader split files. <br>
for the selection analysis, might be better to compare R samples, where leader-containing reads are expected (use HASLEAD) to X samples, where leader should be absent (use NOLEAD)<br>

next might be useful to realign using split inputs:<br>

Align using leader-specific FASTQs: <br>
*.HASLEAD.fastq.gz for no-Xrn1 (R) <br>
*.NOLEAD.fastq.gz for Xrn1 (X)<br>
see **align_loop.sh** which has slight changes to previous loop 

# 06_subcounts


```bash
ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
mkdir -p "$ROOT/work/06_subcounts"
```

```bash
# readcount -> substitution counts (restricted to degenerate region)

mkdir -p $ROOT/work/06_subcounts
python $ROOT/scripts/readcount_to_subcounts.py \
  --readcount_dir $ROOT/work/05_readcount \
  --out_tsv $ROOT/work/06_subcounts/substitution_counts.tsv \
  --min_depth 1000 \
  --restrict \
  --zikv_deg 31-96 \
  --denv_deg 31-94
```

# 07_DESeq2
treat each substitution feature (`pos_ref>alt`) as a “gene” , use for diff. exp. 

determine - what are the “counts” for each feature? <br>
  – alt_count for pos_ref>alt in each sample<br>

`DESeq2` (in R): compare those counts between X vs R across replicates<br>

input: `substitution_counts.tsv` 
outputs: `deseq2_ZIKV_results.tsv`, `deseq2_DENV_results.tsv`
provides:log2FoldChange, lfcSE, pvalue, padj<br>

then, feed `DESeq2` log2FC into the base-pair scoring step (maintain vs disrupt) <br>

### run the script **deseq2_subs.R** in R (local), using substitution_counts.tsv as input


# 07_l2fc-py

```bash
# substitution counts -> l2fc
python $ROOT/scripts/subcounts_to_l2fc.py \
  --in_tsv $ROOT/work/06_subcounts/substitution_counts.tsv \
  --out_tsv $ROOT/work/06_subcounts/substitution_l2fc.tsv \
  --alpha 1e-6
```

# 08_bp-scores

```bash
# base-pair scores for DENV
python $ROOT/scripts/bp_score_from_l2fc.py \
  --l2fc_tsv $ROOT/work/06_subcounts/substitution_l2fc.tsv \
  --out_tsv  $ROOT/work/06_subcounts/denv_bp_scores.tsv \
  --virus DENV \
  --dotbracket '..........(((((.........((.(((((.((.....))))))))))))))......(((((....))))).......' \
  --seq 'UAAAAGAAGUCAGGCCAUCAUAAAUGCCAUAGCUGGAGUAAACUAUGCAGCCUGUAGCUCCACCUGAGAAGGUGUAAAAAA' \
  --start_pos 24 \
  --deg_start 31 --deg_end 94 \
  --min_subs_per_class 1
```
```bash
# base-pair scores for ZIKV
python $ROOT/scripts/bp_score_from_l2fc.py \
  --l2fc_tsv $ROOT/work/06_subcounts/substitution_l2fc.tsv \
  --out_tsv  $ROOT/work/06_subcounts/zikv_bp_scores.tsv \
  --virus ZIKV \
  --dotbracket '..........(((((((....)).((((((.........))))))..))))).......(((((......)))))........' \
  --seq 'AGCUCAUAGUCAGGCCGAGAACGCCAUGGCACGGAAGAAGCCAUGCUGCCUGUGAGCCCCUCAGAGGACACUGAGUCAAAAAA' \
  --start_pos 25 \
  --deg_start 31 --deg_end 96 \
  --min_subs_per_class 1
  ```