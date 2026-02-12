ROOT=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT
TRIM="$ROOT/work/02_trim"
ALIGN="$ROOT/work/04_align_split"
REFDIR="$ROOT/data/refseqs"

mkdir -p "$ALIGN"
cd "$ALIGN"

for SAMPLEFILE in "$TRIM"/*.fastq.gz; do
  FN=$(basename "$SAMPLEFILE")

  # R samples: keep leader-containing
  if [[ "$FN" == *-R-* && "$FN" == *.HASLEAD.fastq.gz ]]; then
    READS="$SAMPLEFILE"

  # X samples: keep leader-absent
  elif [[ "$FN" == *-X-* && "$FN" == *.NOLEAD.fastq.gz ]]; then
    READS="$SAMPLEFILE"
  else
    continue
  fi

  SAMPLE=$(basename "$READS" .fastq.gz)

  # Choose reference index
  if [[ "$SAMPLE" == *-Z* ]]; then
    REF="$REFDIR/zikv_xr2_wt"
  elif [[ "$SAMPLE" == *-D* ]]; then
    REF="$REFDIR/denv_xr2_wt"
  else
    echo "Unknown virus for $SAMPLE — skipping"
    continue
  fi

  echo "Aligning $SAMPLE"

  bowtie2 -p 8 --very-sensitive-local -x "$REF" -U "$READS" \
    2> "${SAMPLE}.bowtie2.log" \
  | samtools view -bS - \
  | samtools sort -o "${SAMPLE}.sorted.bam"

  samtools index "${SAMPLE}.sorted.bam"
done