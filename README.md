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
└── sc_ad_boilerplate/
    ├── metadata.yaml
    ├── samples.tsv
    └── workflow/
        ├── environment.yml
        ├── 01_fetch_and_fastq.sh
        └── 02_pipeline.sh

# 3. To download data
The pipeline supports two modes:

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
	1.	Converts .sra to paired FASTQ files
	2.	Subsamples reads in a paired-safe manner
	3.	Stores them in:
 ```bash
sc_ad_boilerplate/data/fastq_sub/
```
## Pipeline Execution

```bash
# defaults: RATE=0.10, CHEM=10xv3, THREADS=4
bash sc_ad_boilerplate/workflow/02_pipeline.sh
# if you only want to re-run quantification/analysis using existing fastq_sub/:
bash sc_ad_boilerplate/workflow/02_pipeline.sh
```
# 5. Workflow
## Build reference (kb ref)
Purpose: Prepare kallisto|bustools reference for human ; Tools: kb ; Inputs: none (downloads human reference) ; Outputs: workflow/ref/index.idx, t2g.txt, FASTA files

```bash
kb ref -d human \
  -i sc_ad_boilerplate/workflow/ref/index.idx \
  -g sc_ad_boilerplate/workflow/ref/t2g.txt \
  -f1 sc_ad_boilerplate/workflow/ref/cdna.fa \
  -f2 sc_ad_boilerplate/workflow/ref/introns.fa
```
## Step 2 – Quantification (kb count)
Purpose: Quantify transcripts from paired subsampled FASTQs ; Tools: kb (kallisto|bustools) ; Inputs: subsampled FASTQs from data/fastq_sub/ ; Outputs: matrices under workflow/kb_out/<SRR>/counts_unfiltered/
```bash
kb count \
  -i sc_ad_boilerplate/workflow/ref/index.idx \
  -g sc_ad_boilerplate/workflow/ref/t2g.txt \
  -x 10xv3 -t 4 \
  -o sc_ad_boilerplate/workflow/kb_out/<SRR> \
  sc_ad_boilerplate/data/fastq_sub/<SRR>_1.sub.fastq.gz \
  sc_ad_boilerplate/data/fastq_sub/<SRR>_2.sub.fastq.gz
```
## Step 3 – Analysis (Scanpy)
Purpose: merge samples, QC, UMAP, cluster markers, AD vs Control per-cluster ; Tools: Scanpy ; Inputs: matrices from workflow/kb_out/*/counts_unfiltered/ ; Outputs: figures and tables in workflow/scanpy_out/

```bash
# full pipeline (quantification + analysis + read counts)
bash sc_ad_boilerplate/workflow/02_pipeline.sh
```
# 6. Primary outputs:
1.	sc_ad_boilerplate/workflow/scanpy_out/umap_overview.png – UMAP of all cells
2.	sc_ad_boilerplate/workflow/scanpy_out/cell_counts_by_sample_group.csv – cell counts per sample & group
3.	sc_ad_boilerplate/workflow/scanpy_out/markers_per_cluster_wilcoxon.csv – cluster markers
4.	sc_ad_boilerplate/workflow/scanpy_out/DE_AD_vs_Control_per_cluster.csv – DE results per cluster
5.	sc_ad_boilerplate/workflow/scanpy_out/alz_snrna_merged.h5ad – merged AnnData object

Read counts: the pipeline also prints reads per sample based on R1 (zcat | wc -l / 4).
