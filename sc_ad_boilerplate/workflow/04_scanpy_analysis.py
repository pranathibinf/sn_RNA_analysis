import scanpy as sc
import glob
import numpy as np
import pandas as pd
import anndata as ad

# --- Load all kb count outputs ---
sample_dirs = glob.glob("sc_ad_boilerplate/workflow/kb_out/*/counts_unfiltered")
adatas = []
for d in sample_dirs:
    adata_tmp = sc.read_mtx(f"{d}/spliced.mtx")
    adata_tmp.var_names = [x.strip() for x in open(f"{d}/spliced.genes.txt")]
    adata_tmp.obs_names = [x.strip() for x in open(f"{d}/spliced.barcodes.txt")]
    adata_tmp = adata_tmp.T  # transpose → cells × genes
    sample_id = d.split("/")[-2]
    adata_tmp.obs["sample"] = sample_id
    adatas.append(adata_tmp)

# --- Merge into single AnnData ---
adata = ad.concat(adatas, join="outer", label="sample_id", keys=[a.obs["sample"][0] for a in adatas])

# --- QC ---
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_cells(adata, max_counts=50000)
adata = adata[adata.obs.pct_counts_mt < 15]

# --- Normalize, HVGs ---
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.scale(adata, max_value=10)

# --- PCA, neighbors, UMAP ---
n_comps = min(30, adata.n_vars, adata.n_obs) - 1
n_comps = max(2, n_comps)   # ensure at least 2 PCs
print(f"Using n_comps={n_comps} for PCA")

sc.tl.pca(adata, svd_solver="arpack", n_comps=n_comps)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=n_comps)
sc.tl.umap(adata)

# --- Clustering ---
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)
import os

# --- Ensure output directory exists ---
out_dir = "sc_ad_boilerplate/workflow/scanpy_out"
os.makedirs(out_dir, exist_ok=True)

# --- Save outputs ---
sc.pl.umap(adata, color=["sample", "leiden"], wspace=0.4, save="_overview.png")
adata.write(f"{out_dir}/alz_snrna_merged.h5ad")
