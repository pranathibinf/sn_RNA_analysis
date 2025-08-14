!/usr/bin/env bash
# Stream SRR → paired FASTQ.gz with subsampling (paired-safe)
# Default: streaming (no .sra on disk). Optionally USE_PREFETCH=1 for local .sra.
# Vars: RATE (default 0.10), THREADS (default 4), USE_PREFETCH (default 0)
set -euo pipefail

RATE="${RATE:-0.10}"
THREADS="${THREADS:-4}"
USE_PREFETCH="${USE_PREFETCH:-0}"

ROOT="sc_ad_boilerplate"
FASTQ_DIR="$ROOT/data/fastq_sub"
SRA_DIR="$ROOT/data/sra"
SAMPLES_TSV="$ROOT/samples.tsv"

mkdir -p "$FASTQ_DIR"
[[ "$USE_PREFETCH" -eq 1 ]] && mkdir -p "$SRA_DIR"

# checks
for x in fastq-dump reformat.sh; do
  command -v "$x" >/dev/null || { echo "Missing: $x"; exit 1; }
done
[[ "$USE_PREFETCH" -eq 1 ]] && command -v prefetch >/dev/null || true

echo "[stream] Subsample rate: $RATE"
tail -n +2 "$SAMPLES_TSV" | while read -r srr group; do
  [[ -z "$srr" ]] && continue
  out1="$FASTQ_DIR/${srr}_1.sub.fastq.gz"
  out2="$FASTQ_DIR/${srr}_2.sub.fastq.gz"
  if [[ -s "$out1" && -s "$out2" ]]; then
    echo "  - $srr: FASTQs exist, skip"
    continue
  fi

  if [[ "$USE_PREFETCH" -eq 1 ]]; then
    echo "  - $srr: prefetch → local .sra → convert"
    command -v prefetch >/dev/null || { echo "Missing: prefetch"; exit 1; }
    prefetch "$srr" -O "$SRA_DIR"
    fastq-dump --stdout --split-spot --skip-technical "$SRA_DIR/$srr/$srr.sra" \
      | reformat.sh in=stdin.fq int=t samplerate="$RATE" overwrite=t \
          out1="$out1" out2="$out2" gzip=9
  else
    echo "  - $srr: streaming directly"
    fastq-dump --stdout --split-spot --skip-technical "$srr" \
      | reformat.sh in=stdin.fq int=t samplerate="$RATE" overwrite=t \
          out1="$out1" out2="$out2" gzip=9
  fi
done

echo "[done] FASTQs written to $FASTQ_DIR"
