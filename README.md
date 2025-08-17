# 1. Data sources
This task is based on publicly available sequencing data from a study of Alzheimer’s Disease and Down Syndrome (https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA975472) .The dataset includes multiple samples under different conditions (AD vs Control) and was originally sequenced using Illumina paired-end 2×150. The subsampled FASTQs are stored in sc_ad_boilerplate/data/ upon running the pipeline and are used as the inputs for the workflow.

	• BioProject: PRJNA975472
	• GEO SuperSeries: GSE233208
	• Paper: Miyoshi et al., Nature Genetics (2024)
    • Assay: single-nucleus RNA-seq (10x)

## Subset used:
| SRR        | Group    |
|-----------:|----------|
| SRR24710554 | Control |
| SRR24710556 | AD      |
| SRR24710558 | Control |
| SRR24710560 | AD      | 

## 2. Repository layout:

sc-alzheimer-analysis/
├── README.md
├── LICENSE
└── sc_ad_boilerplate/
    ├── metadata.yaml
    ├── samples.tsv
    └── workflow/
        ├── 01_fetch_and_fastq.sh
        ├── 02_build_ref.sh
        ├── 03_kb_count.sh
        ├── 04_scanpy_analysis.py
        └── environment.yml

# 3. Installation
Create and activate the environment:

- **Prefetch**: download `.sra` locally, then convert to FASTQ  
- **Stream**: download and process directly without storing `.sra` locally

## Create and activate environment
```bash
conda env create -f sc_ad_boilerplate/workflow/environment.yml
conda activate ad_snrna

# Download and convert to paired, subsampled FASTQ.gz (default: 10% subsample)
bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh
```
```bash
# stream directly from SRA without saving .sra locally
STREAM=1 bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh
# change subsample rate (e.g., 20%)
RATE=0.20 bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh
```

# 4.  Pre-processing / subsampling
The script:
  1. Converts .sra to paired FASTQ files
  2. Subsamples reads in a paired-safe manner (default: 10%)
  3. Stores them in:
 ```bash
sc_ad_boilerplate/data/fastq_sub/
```
## Pipeline Execution

```bash
# defaults: RATE=0.10, THREADS=4
bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh

# if you only want to stream directly from SRA without storing .sra:
STREAM=1 bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh

# change subsample rate (e.g., 20%)
RATE=0.20 bash sc_ad_boilerplate/workflow/01_fetch_and_fastq.sh
```
# 5. Workflow
## Step 1 – Build reference (kb ref) 
Purpose: Prepare kallisto|bustools reference for **mouse (Mus musculus, GRCm38, Ensembl 98)**  
Tools: kb  ;
Inputs: Ensembl FASTA + GTF (downloaded automatically)  ;
Outputs: workflow/ref/index.idx, t2g.txt, FASTA files

```bash
bash sc_ad_boilerplate/workflow/02_build_ref.sh
```
## Step 2 – Quantification (kb count)
Purpose: Quantify transcripts from paired subsampled FASTQs  ;
Tools: kb (kallisto|bustools)  ;
Inputs: subsampled FASTQs from `data/fastq_sub/`  ;
Outputs: matrices under `workflow/kb_out/<SRR>/counts_unfiltered/`workflow/kb_out/<SRR>/counts_unfiltered/
```bash
bash sc_ad_boilerplate/workflow/03_kb_count.sh
```
## Step 3 – Analysis (Scanpy)
Purpose: Merge samples, perform QC, UMAP, clustering, and compare AD vs Control  ;
Tools: Scanpy (v1.10.2)  ;
Inputs: matrices from `workflow/kb_out/*/counts_unfiltered/`  ;
Outputs: figures and tables in `workflow/scanpy_out/`

```bash
python sc_ad_boilerplate/workflow/04_scanpy_analysis.py
```
# 6. Primary outputs:
1. sc_ad_boilerplate/workflow/scanpy_out/umap_overview.png – UMAP of all cells
2. sc_ad_boilerplate/workflow/scanpy_out/cell_counts_by_sample_group.csv – cell counts per sample & group
3. sc_ad_boilerplate/workflow/scanpy_out/markers_per_cluster_wilcoxon.csv – cluster markers
4. sc_ad_boilerplate/workflow/scanpy_out/DE_AD_vs_Control_per_cluster.csv – DE results per cluster
5. sc_ad_boilerplate/workflow/scanpy_out/alz_snrna_merged.h5ad – merged AnnData object

**Read counts check:**
```bash
zcat sc_ad_boilerplate/data/fastq_sub/<SRR>_1.sub.fastq.gz | wc -l
```
