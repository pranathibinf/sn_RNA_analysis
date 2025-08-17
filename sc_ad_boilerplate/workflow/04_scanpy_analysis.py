import scanpy as sc
import glob
import numpy as np
import pandas as pd

# Load all kb count outputs
sample_dirs = glob.glob("sc_ad_boilerplate/workflow/kb_out/*/counts_unfiltered")
adatas = []
for d in sample_dirs:
    ad = sc.read_mtx(f"{d}/spliced.mtx")
    ad.var_names = [x.strip() for x in open(f"{d}/spliced.genes.txt")]
    ad.obs_names = [x.strip() for x in open(f"{d}/spliced.barcodes.txt")]
    ad = ad.T
    sample_id = d.split("/")[-2]
    ad.obs["sample"] = sample_id
    adatas.append(ad)

# Merge into single AnnData
adata = adatas[0].concatenate(
    adatas[1:], join="outer", batch_key="sample_id",
    batch_categories=[a.obs["sample"][0] for a in adatas]
)

# QC
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_cells(adata, max_counts=50000)
adata = adata[adata.obs.pct_counts_mt < 15]

# Normalize, HVGs, PCA, UMAP
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=True)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=30)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=1.0, flavor="igraph", n_iterations=2, directed=False)

# Save outputs
sc.pl.umap(adata, color=["sample", "leiden"], wspace=0.4, save="_overview.png")
adata.write("sc_ad_boilerplate/workflow/scanpy_out/alz_snrna_merged.h5ad")
