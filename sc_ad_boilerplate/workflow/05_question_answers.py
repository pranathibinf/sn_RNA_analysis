# --- Imports ---
import scanpy as sc
import numpy as np
import glob, os

# ========= Load all samples =========
base_dir = "sc_ad_boilerplate/workflow/kb_out"
sample_dirs = glob.glob(f"{base_dir}/*/counts_unfiltered")

print("Found sample dirs:", sample_dirs)
if len(sample_dirs) == 0:
    raise FileNotFoundError("No kb_out/counts_unfiltered folders found. Check your path!")

adatas = []
for d in sample_dirs:
    print(f"Loading {d}")
    ad = sc.read_mtx(f"{d}/spliced.mtx")
    ad.var_names = [x.strip() for x in open(f"{d}/spliced.genes.txt")]
    ad.obs_names = [x.strip() for x in open(f"{d}/spliced.barcodes.txt")]
    ad = ad.T
    sample_id = d.split("/")[-2]
    ad.obs["sample"] = sample_id
    adatas.append(ad)

# Merge AnnData objects
if len(adatas) == 1:
    adata = adatas[0]
    print("Only one sample detected — using it directly.")
else:
    adata = sc.concat(adatas, join="outer", label="sample_id", keys=[a.obs["sample"][0] for a in adatas])
    print(f"Merged {len(adatas)} samples")

# ========= Q1: gene with highest counts =========
gene_counts = np.array(adata.X.sum(axis=0)).flatten()
top_gene = adata.var_names[np.argmax(gene_counts)]
print("Q1:", top_gene)

# ========= Q2: highly variable genes =========
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat", subset=False)
hvg_count = int(adata.var["highly_variable"].sum())
print("Q2:", hvg_count)

# ========= Q3 + Q4: PCA variance =========
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack")

pc_var = adata.uns["pca"]["variance_ratio"]
pc1_var = pc_var[0] * 100
print("Q3:", round(pc1_var, 2), "%")

# Rank top 3 PCs
top3 = np.argsort(-pc_var)[:3] + 1
print("Q4:", [f"PC{pc}" for pc in top3])

# ========= Q5: clustering =========
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=10)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)

cluster_sizes = adata.obs["leiden"].value_counts()
largest_cluster = cluster_sizes.idxmax()
print("Q5:", largest_cluster)

# ========= Save outputs =========
os.makedirs("sc_ad_boilerplate/workflow/scanpy_out", exist_ok=True)
adata.write("sc_ad_boilerplate/workflow/scanpy_out/alz_snrna_merged.h5ad")
sc.pl.umap(adata, color=["sample", "leiden"], wspace=0.4, save="_overview.png")
