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
