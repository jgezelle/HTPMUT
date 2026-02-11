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
### (from `RAW` dir):
```bash
ln -s $SRC/A3413_R1.fastq.gz  260120-r0400-DR-01-A3413_R1.fastq.gz
ln -s $SRC/A3429_R1.fastq.gz  260120-r0400-DX-01-A3429_R1.fastq.gz

ln -s $SRC/A3419_R1.fastq.gz  260120-r0400-DR-02-A3419_R1.fastq.gz
ln -s $SRC/A3499_R1.fastq.gz  260120-r0400-DX-02-A3499_R1.fastq.gz

ln -s $SRC/A3485_R1.fastq.gz  260120-r0400-DR-03-A3485_R1.fastq.gz
ln -s $SRC/A3484_R1.fastq.gz  260120-r0400-DX-03-A3484_R1.fastq.gz

ln -s $SRC/A3483_R1.fastq.gz  260120-r0399-ZR-01-A3483_R1.fastq.gz
ln -s $SRC/A3466_R1.fastq.gz  260120-r0399-ZX-01-A3466_R1.fastq.gz

ln -s $SRC/A3456_R1.fastq.gz  260120-r0399-ZR-02-A3456_R1.fastq.gz
ln -s $SRC/A3452_R1.fastq.gz  260120-r0399-ZX-02-A3452_R1.fastq.gz

ln -s $SRC/A3451_R1.fastq.gz  260120-r0399-ZR-03-A3451_R1.fastq.gz
ln -s $SRC/A3450_R1.fastq.gz  260120-r0399-ZX-03-A3450_R1.fastq.gz

```

**check symlinks worked:**<br>
```bash
ls -lh
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
conda install -c bioconda bedtools
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
### try pileup without read depth limitation (takes more memory but preserves depth): EDIT: memory hog, doesn't work well
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
