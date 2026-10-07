# Automated End-to-End RNA-Seq Pipeline (Transcriptomics)

An automated, production-ready transcriptomics pipeline that processes raw RNA-Seq data (`.fastq`/`.fastq.gz`) into publication-ready differential expression data and pathway enrichment visualizations. 

This workflow bridges standard Linux command-line tools for sequence alignment with R Bioconductor packages for statistical modeling, handling both paired-end and single-end reads automatically.

## Pipeline Architecture
The workflow is divided into two phases: a Linux/Bash pre-processing phase and an R-based statistical analysis phase.

**1. Data Processing (Bash)**
*   **Quality Control & Trimming:** `fastp` (Handles both compressed and uncompressed FASTQ files)
*   **Splice-Aware Alignment:** `HISAT2` (Auto-builds indices from reference `.fa` genomes)
*   **Read Quantification:** `featureCounts` (Generates standard count matrices against `.gtf` annotations)
*   **Data Wrangling:** `samtools` & `sed` (Automated BAM sorting and regex-based matrix cleanup)

**2. Statistical Modeling (R / Bioconductor)**
*   **Library Normalization & Dispersion:** `DESeq2` (Models data via Negative Binomial distribution)
*   **Differential Expression:** `EnhancedVolcano` & `pheatmap`
*   **Gene Ontology (GO) & Pathway Enrichment:** `clusterProfiler` & `msigdbr`
*   **Available Modules:** Domain-specific scripts provided for Human (`org.Hs.eg.db`), Bacteria (`org.EcK12.eg.db`), and Plant (`org.At.tair.db`) transcriptomes.

## Key Engineering Features
*   **Idempotent Execution:** Bash scripts are configured to prevent file-overwriting errors if re-run, checking local directories before pulling master references.
*   **Fail-Safes & Error Trapping:** R scripts include matrix-rank checks to prevent silent statistical failures if biological conditions lack variation, and dataframes are forced to bypass base R's hyphen-conversion (check.names=FALSE) to ensure count matrices sync perfectly with metadata labels.
*   **Isolated Environment:** Fully reproducible via a unified Conda `environment.yml`.

## How to Run the Pipeline

### Setup the Environment
Clone this repository and build the dedicated Conda environment to resolve all dependencies automatically:
```bash
conda env create -f environment_rnaseq.yml
conda activate rnaseq_env
