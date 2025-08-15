#!/usr/bin/env bash
set -euo pipefail

# Run from repo root (this file lives in workflow/)
cd "$(dirname "$0")/.."

# -------- config (overridable) --------
FASTQ_DIR="${FASTQ_DIR:-sc_ad_boilerplate/data/fastq_sub}"   # input FASTQs
SAMPLES_TSV="${SAMPLES_TSV:-workflow/samples.tsv}"           # SRR/group table
REF_DIR="workflow/ref"
KB_OUT_DIR="workflow/kb_out"
SCANPY_OUT_DIR="workflow/scanpy_out"
THREADS="${THREADS:-4}"
CHEM="${CHEM:-10xv3}"
GZCAT_CMD="gunzip -c"

log(){ printf "[%s] %s\n" "$(date '+%H:%M:%S')" "$*"; }
die(){ echo "ERROR: $*" >&2; exit 1; }
need(){ command -v "$1" >/dev/null 2>&1 || die "Missing '$1' in PATH. Activate env?"; }

need kb
need python3
mkdir -p "$REF_DIR" "$KB_OUT_DIR" "$SCANPY_OUT_DIR" "$FASTQ_DIR"

# -------- samples.tsv  --------
if [[ ! -f "$SAMPLES_TSV" ]]; then
  log "samples.tsv not found — creating one from FASTQ names"
  SRRS=$(ls "$FASTQ_DIR"/*_1.sub.fastq.gz 2>/dev/null | sed -E 's#.*/(SRR[0-9]+)_1\.sub\.fastq\.gz#\1#' | sort -u || true)
  if [[ -z "${SRRS:-}" ]]; then
    SRRS="SRR24710554 SRR24710556 SRR24710558 SRR24710560"
  fi
  {
    echo -e "SRR\tgroup"
    for s in $SRRS; do
      case "$s" in
        SRR24710554|SRR24710558) grp="Control" ;;
        SRR24710556|SRR24710560) grp="AD" ;;
        *) grp="Unknown" ;;
      esac
      echo -e "$s\t$grp"
    done
  } > "$SAMPLES_TSV"
  log "wrote $SAMPLES_TSV"
fi

# -------- reference  --------
REF_IDX="$REF_DIR/index.idx"
if [[ ! -f "$REF_IDX" ]]; then
  log "ref: building kallisto|bustools human reference"
  set +e
  kb ref -d human \
    -i "$REF_IDX" \
    -g "$REF_DIR/t2g.txt" \
    -f1 "$REF_DIR/cdna.fa" \
    -f2 "$REF_DIR/introns.fa"
  rc=$?
  set -e
  if [[ $rc -ne 0 || ! -f "$REF_IDX" ]]; then
    log "ref: primary failed; using fallback download"
    URL="https://github.com/pachterlab/kallisto-transcriptome-indices/releases/download/v1/human_index_standard.tar.xz"
    TGT="workflow/tmp_human_index_standard.tar.xz"
    mkdir -p workflow
    curl -L --fail --retry 3 --retry-delay 3 -o "$TGT" "$URL"
    tar -xJf "$TGT" -C "$REF_DIR"
  fi
  [[ -f "$REF_IDX" ]] || die "reference build failed"
else
  log "ref: found $REF_IDX"
fi

# -------- quantification --------
log "kb: count (chem=$CHEM, threads=$THREADS)"
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r srr group; do
  [[ "$srr" == SRR* ]] || continue
  fq1="$FASTQ_DIR/${srr}_1.sub.fastq.gz"
  fq2="$FASTQ_DIR/${srr}_2.sub.fastq.gz"
  out="$KB_OUT_DIR/$srr"
  if [[ -s "$fq1" && -s "$fq2" ]]; then
    log "kb: $srr"
    kb count \
      -i "$REF_IDX" \
      -g "$REF_DIR/t2g.txt" \
      -x "$CHEM" \
      -t "$THREADS" \
      -o "$out" \
      "$fq1" "$fq2"
  else
    log "warn: missing FASTQs for $srr in $FASTQ_DIR"
  fi
done

# -------- 10x-compat shim (gz output) --------
log "compat: creating 10x-style files"
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r srr group; do
  [[ "$srr" == SRR* ]] || continue
  d="$KB_OUT_DIR/$srr/counts_unfiltered"
  [[ -d "$d" ]] || { log "compat: skip $srr (no counts_unfiltered)"; continue; }

  # matrix.mtx.gz
  if [[ -f "$d/cells_x_genes.mtx" ]]; then
    gzip -c "$d/cells_x_genes.mtx" > "$d/matrix.mtx.gz"
  fi

  # barcodes.tsv + barcodes.tsv.gz
  if [[ -f "$d/cells_x_genes.barcodes.txt" ]]; then
    ln -sf "cells_x_genes.barcodes.txt" "$d/barcodes.tsv"
    gzip -c "$d/cells_x_genes.barcodes.txt" > "$d/barcodes.tsv.gz"
  fi

  # features.tsv.gz (gene_id \t gene_symbol \t Gene Expression)
  if [[ -f "$d/cells_x_genes.genes.txt" ]]; then
    if [[ -f "$d/cells_x_genes.genes.names.txt" ]]; then
      paste "$d/cells_x_genes.genes.txt" "$d/cells_x_genes.genes.names.txt" \
        | awk -F'\t' 'BEGIN{OFS="\t"}{sym=$2; if(sym==""){sym=$1}; print $1, sym, "Gene Expression"}' \
        | gzip -c > "$d/features.tsv.gz"
    else
      awk '{print $1"\t"$1"\tGene Expression"}' "$d/cells_x_genes.genes.txt" \
        | gzip -c > "$d/features.tsv.gz"
    fi
  fi
done

# -------- Scanpy analysis  --------
log "scanpy: analysis"
python3 <<'PYCODE'
import os, gzip, numpy as np, pandas as pd
from scipy.io import mmread
from scipy.sparse import csr_matrix
import anndata as ad
import scanpy as sc

KB_OUT_DIR = "workflow/kb_out"
SCANPY_OUT_DIR = "workflow/scanpy_out"
SAMPLES_TSV = "workflow/samples.tsv"

os.makedirs(SCANPY_OUT_DIR, exist_ok=True)
sc.settings.figdir = SCANPY_OUT_DIR

def read_lines(path):
    if path.endswith(".gz"):
        return [l.decode().strip() for l in gzip.open(path, "rb")]
    return [l.strip() for l in open(path, "r")]

samples = pd.read_csv(SAMPLES_TSV, sep="\t")
adatas = []
for _, row in samples.iloc[1:].iterrows() if list(samples.columns)!=['SRR','group'] else samples.iloc[1:].iterrows():
    srr = str(row["SRR"])
    d = os.path.join(KB_OUT_DIR, srr, "counts_unfiltered")

    mtx = os.path.join(d, "matrix.mtx.gz")
    if not os.path.exists(mtx): mtx = os.path.join(d, "matrix.mtx")
    feats = os.path.join(d, "features.tsv.gz")
    if not os.path.exists(feats): feats = os.path.join(d, "features.tsv")
    bar = os.path.join(d, "barcodes.tsv.gz")
    if not os.path.exists(bar): bar = os.path.join(d, "barcodes.tsv")

    if not (os.path.exists(mtx) and os.path.exists(feats) and os.path.exists(bar)):
        print(f"[scanpy] skip {srr}: missing matrix/features/barcodes in {d}")
        continue

    X = mmread(mtx).tocsr()
    features = [ln.split("\t") for ln in read_lines(feats)]
    gene_ids = [f[0] for f in features]
    gene_syms = [ (f[1] if len(f)>1 and f[1] else f[0]) for f in features ]
    barcodes = read_lines(bar)

    # Orient matrix: cells x genes
    if X.shape[0] == len(gene_ids) and X.shape[1] == len(barcodes):
        X = X.transpose().tocsr()
    elif X.shape[0] == len(barcodes) and X.shape[1] == len(gene_ids):
        pass
    else:
        raise RuntimeError(f"{srr}: matrix shape {X.shape} doesn't match genes={len(gene_ids)} / cells={len(barcodes)}")

    obs = pd.DataFrame(index=barcodes)
    obs["barcode"] = barcodes
    obs["sample"] = srr
    obs["group"] = row.get("group", "Unknown")

    var = pd.DataFrame(index=gene_ids)
    var["gene_id"] = gene_ids
    var["gene_symbol"] = gene_syms
    var_names = pd.Series(gene_syms).replace("", np.nan).fillna(pd.Series(gene_ids)).values

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names = var_names
    adatas.append(adata)

if not adatas:
    raise SystemExit("[scanpy] No matrices loaded — nothing to analyze.")

adata = ad.concat(adatas, join="outer", label="sample_concat", keys=[a.obs["sample"][0] for a in adatas])

# Basic pipeline
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
sc.tl.umap(adata)

# Try Leiden; if igraph/leidenalg missing, skip if not
leiden_ok = True
try:
    import igraph, leidenalg  # noqa: F401
    sc.tl.leiden(adata, key_added="leiden")
except Exception as e:
    print(f"[scanpy] Leiden skipped (install python-igraph & leidenalg to enable). Reason: {e}")
    leiden_ok = False

# Plots & save
colors = ['group']
if leiden_ok: colors.append('leiden')
sc.pl.umap(adata, color=colors, save="_overview.png", show=False)
adata.write(os.path.join(SCANPY_OUT_DIR, "alz_snrna_merged.h5ad"))
print(f"[scanpy] wrote {os.path.join(SCANPY_OUT_DIR, 'alz_snrna_merged.h5ad')}")
PYCODE

# -------- read-count summary (R1) --------
printf "sample\treads\n"
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r srr group; do
  [[ "$srr" == SRR* ]] || continue
  fq1="$FASTQ_DIR/${srr}_1.sub.fastq.gz"
  [[ -s "$fq1" ]] || continue
  reads=$($GZCAT_CMD "$fq1" | wc -l | awk '{printf "%.0f", $1/4}')
  printf "%s\t%s\n" "$srr" "$reads"
done

log "done: pipeline completed."

