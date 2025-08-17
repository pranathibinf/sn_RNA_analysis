#!/bin/bash
set -euo pipefail

THREADS=4
REF_DIR="sc_ad_boilerplate/workflow/ref"
FASTQ_DIR="sc_ad_boilerplate/data/fastq_sub"
OUT_DIR="sc_ad_boilerplate/workflow/kb_out"

mkdir -p "$OUT_DIR"

SAMPLES=($(ls $FASTQ_DIR/*_1.sub.fastq.gz | sed 's/.*\///; s/_1\.sub\.fastq\.gz//' | sort -u))

for SAMPLE in "${SAMPLES[@]}"; do
  echo "[INFO] Processing $SAMPLE"

  kb count \
    -i "$REF_DIR/index.idx" \
    -g "$REF_DIR/t2g.txt" \
    -x 10xv3 \
    -t "$THREADS" \
    -o "$OUT_DIR/$SAMPLE" \
    --workflow lamanno \
    -c1 "$REF_DIR/cdna_t2c.txt" \
    -c2 "$REF_DIR/intron_t2c.txt" \
    "$FASTQ_DIR/${SAMPLE}_1.sub.fastq.gz" \
    "$FASTQ_DIR/${SAMPLE}_2.sub.fastq.gz"

  echo "[INFO] Completed $SAMPLE"
  echo "----------------------------------------"
done
