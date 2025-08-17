#!/bin/bash
set -euo pipefail

THREADS=4
RATE=${RATE:-0.10}
FASTQ_DIR="sc_ad_boilerplate/data/fastq_sub"
SRA_DIR="sc_ad_boilerplate/data/sra"

mkdir -p "$FASTQ_DIR" "$SRA_DIR"

# SRRs (Control + AD)
SAMPLES=(
  SRR24710554
  SRR24710556
  SRR24710558
  SRR24710560
)

i=0; n=${#SAMPLES[@]}
for acc in "${SAMPLES[@]}"; do
  i=$((i+1))
  echo "[$i/$n] Processing $acc"

  # Prefetch SRA
  if [[ "${STREAM:-0}" -eq 0 ]]; then
    prefetch -O "$SRA_DIR" "$acc"
    sra="$SRA_DIR/${acc}.sra"
  else
    sra="$acc"
  fi

  # Convert to FASTQ + subsample
  fastq-dump --stdout --split-spot --skip-technical "$sra" \
  | seqtk sample -s100 - "$RATE" \
  | paste - - - - - - - - \
  | tee >(awk '{print $1 ORS $2 ORS $3 ORS $4}' | tr '\t' '\n' \
          | gzip > "${FASTQ_DIR}/${acc}_1.sub.fastq.gz") \
        >(awk '{print $5 ORS $6 ORS $7 ORS $8}' | tr '\t' '\n' \
          | gzip > "${FASTQ_DIR}/${acc}_2.sub.fastq.gz") \
  > /dev/null
done
