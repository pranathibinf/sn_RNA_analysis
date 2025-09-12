#!/bin/bash
set -euo pipefail

REF_DIR="sc_ad_boilerplate/workflow/ref"
mkdir -p "$REF_DIR"

# Mouse GRCm38, Ensembl 98
wget -nc ftp://ftp.ensembl.org/pub/release-98/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -P "$REF_DIR"
wget -nc ftp://ftp.ensembl.org/pub/release-98/gtf/mus_musculus/Mus_musculus.GRCm38.98.gtf.gz -P "$REF_DIR"

kb ref \
  -i "$REF_DIR/index.idx" \
  -g "$REF_DIR/t2g.txt" \
  -f1 "$REF_DIR/cdna.fa" \
  -f2 "$REF_DIR/intron.fa" \
  -c1 "$REF_DIR/cdna_t2c.txt" \
  -c2 "$REF_DIR/intron_t2c.txt" \
  --workflow lamanno \
  "$REF_DIR/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz" \
  "$REF_DIR/Mus_musculus.GRCm38.98.gtf.gz"
