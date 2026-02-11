# Directory with trimmed and split files
TRIM=/groups/as6282_gp/scratch_bkup/jgg2144/HTPMUT/work/02_trim

# Output summary table
SUMMARY="$TRIM/leader_summary.tsv"

# Header
echo -e "Sample\tHASLEAD\tNOLEAD\tTotal\tLeaderFraction" > "$SUMMARY"

# Loop over all HASLEAD files
for HAS in "$TRIM"/*.HASLEAD.fastq.gz; do
  # Base sample name
  BASE=$(basename "$HAS" .HASLEAD.fastq.gz)
  
  # Corresponding NOLEAD file
  NOLEAD="$TRIM/${BASE}.NOLEAD.fastq.gz"

  # Count reads by line count /4
  HAS_COUNT=$(( $(zcat "$HAS" | wc -l) / 4 ))
  NO_COUNT=$(( $(zcat "$NOLEAD" | wc -l) / 4 ))

  # Total reads
  TOTAL=$(( HAS_COUNT + NO_COUNT ))

  # Fraction of reads retaining leader (decimal)
  # bc used for floating point division
  FRACTION=$(echo "scale=4; $HAS_COUNT / $TOTAL" | bc)

  # Append to summary table
  echo -e "${BASE}\t${HAS_COUNT}\t${NO_COUNT}\t${TOTAL}\t${FRACTION}" >> "$SUMMARY"
done

echo "Summary table written to $SUMMARY"
