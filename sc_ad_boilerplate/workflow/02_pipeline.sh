#!/usr/bin/env bash
# Quantification (kb) + Scanpy analysis + read counts
# Vars: THREADS (default 4), CHEM (default 10xv3)
set -euo pipefail

THREADS="${THREADS:-4}"
CHEM="${CHEM:-10xv3}"

ROOT="sc_ad_boilerplate"
FASTQ_DIR="$ROOT/data/fastq_sub"
REF_DIR="$ROOT/workflow/ref"
KB_OUT="$ROOT/workflow/kb_out"
SCANPY_OUT="$ROOT/workflow/scanpy_out"
SAMPLES_TSV="$ROOT/samples.tsv"

mkdir -p "$REF_DIR" "$KB_OUT" "$SCANPY_OUT"

for x in kb python; do
  command -v "$x" >/dev/null || { echo "Missing: $x"; exit 1; }
done

# build reference once
if [[ ! -s "$REF_DIR/index.idx" ]]; then
  echo "[ref] kb ref (human)"
  kb ref -d human -i "$REF_DIR/index.idx" -g "$REF_DIR/t2g.txt" \
        -f1 "$REF_DIR/cdna.fa" -f2 "$REF_DIR/introns.fa"
else
  echo "[ref] found $REF_DIR/index.idx"
fi

# kb count per sample
echo "[kb] count (chem=$CHEM, threads=$THREADS)"
tail -n +2 "$SAMPLES_TSV" | while read -r srr group; do
  [[ -z "$srr" ]] && continue
  r1="$FASTQ_DIR/${srr}_1.sub.fastq.gz"
  r2="$FASTQ_DIR/${srr}_2.sub.fastq.gz"
  out="$KB_OUT/$srr"
  if [[ ! -s "$r1" || ! -s "$r2" ]]; then
    echo "[warn] missing FASTQs for $srr"; continue
  fi
  if [[ -d "$out/counts_unfiltered" ]]; then
    echo "  - $srr: exists"; continue
  fi
  kb count -i "$REF_DIR/index.idx" -g "$REF_DIR/t2g.txt" -x "$CHEM" -t "$THREADS" -o "$out" "$r1" "$r2"
done

# read counts (R1)
echo -e "sample\treads"
for r1 in "$FASTQ_DIR"/*_1.sub.fastq.gz; do
  s=$(basename "$r1" _1.sub.fastq.gz)
  n=$(zcat "$r1" | wc -l | awk '{print int($1/4)}')
  echo -e "$s\t$n"
done

# Scanpy: merge, QC, UMAP, markers, DE per cluster
python - <<'PY'
import os, pandas as pd, scanpy as sc
from scipy.io import mmread

ROOT="sc_ad_boilerplate"
SAMPLES=os.path.join(ROOT,"samples.tsv")
KB=os.path.join(ROOT,"workflow","kb_out")
OUT=os.path.join(ROOT,"workflow","scanpy_out")
os.makedirs(OUT, exist_ok=True)

df=pd.read_csv(SAMPLES, sep=None, engine="python")

def read_one(srr):
    d=os.path.join(KB,srr,"counts_unfiltered")
    m=os.path.join(d,"cells_x_genes.mtx")
    b=os.path.join(d,"cells_x_genes.barcodes.txt")
    g1=os.path.join(d,"cells_x_genes.genes.names.txt")
    g2=os.path.join(d,"cells_x_genes.genes.txt")
    if not os.path.exists(m): return None
    X=mmread(m).tocsr()
    ad=sc.AnnData(X)
    ad.obs_names=pd.read_csv(b,header=None)[0].astype(str).values
    ad.var_names=pd.read_csv(g1 if os.path.exists(g1) else g2,header=None)[0].astype(str).values
    return ad

ads=[]
for _,row in df.iterrows():
    srr=row['srr']; grp=row['group']
    ad=read_one(srr)
    if ad is None: 
        print(f"[WARN] missing counts for {srr}")
        continue
    ad.obs['sample']=srr; ad.obs['group']=grp; ad.obs['batch']=srr
    ads.append(ad)

assert ads, "No matrices found"

adata=ads[0].concatenate(ads[1:], batch_key="concat_batch",
                         batch_categories=[a.obs['sample'][0] for a in ads],
                         index_unique="-")
adata.var["mt"]=adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)

adata=adata[adata.obs["n_genes_by_counts"]>=100].copy()
adata=adata[adata.obs["n_genes_by_counts"]<=8000].copy()
adata=adata[adata.obs["pct_counts_mt"]<25].copy()
sc.pp.filter_genes(adata, min_cells=3)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.tl.pca(adata, svd_solver="arpack")
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

sc.pl.umap(adata, color=["group","sample","leiden","pct_counts_mt"], wspace=0.4, show=False)
import matplotlib.pyplot as plt
plt.savefig(os.path.join(OUT,"umap_overview.png"), dpi=200, bbox_inches="tight")

adata.obs.groupby(["sample","group"]).size().to_frame("n_cells").reset_index()\
  .to_csv(os.path.join(OUT,"cell_counts_by_sample_group.csv"), index=False)

sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon")
sc.get.rank_genes_groups_df(adata, group=None)\
  .to_csv(os.path.join(OUT,"markers_per_cluster_wilcoxon.csv"), index=False)

outs=[]
for cl in sorted(adata.obs["leiden"].unique(), key=str):
    ad_cl=adata[adata.obs["leiden"]==cl].copy()
    c=ad_cl.obs["group"].value_counts().to_dict()
    if c.get("AD",0)>=5 and c.get("Control",0)>=5:
        sc.tl.rank_genes_groups(ad_cl,"group",groups=["AD"],reference="Control",method="wilcoxon")
        df=sc.get.rank_genes_groups_df(ad_cl, group="AD")
        df.insert(0,"cluster",cl); df.insert(1,"n_AD",c.get("AD",0)); df.insert(2,"n_Control",c.get("Control",0))
        outs.append(df)
    else:
        outs.append(pd.DataFrame({"cluster":[cl],"n_AD":[c.get("AD",0)],"n_Control":[c.get("Control",0)]}))
if outs:
    pd.concat(outs, ignore_index=True).to_csv(os.path.join(OUT,"DE_AD_vs_Control_per_cluster.csv"), index=False)

adata.write(os.path.join(OUT,"alz_snrna_merged.h5ad"))
print("[scanpy] results in", OUT)
PY

echo "[ok] Pipeline complete. See $SCANPY_OUT"
