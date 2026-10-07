#!/bin/bash
set -e

echo "====================================================================="
echo " UNIVERSAL PAIRED-END RNA-SEQ PIPELINE"
echo "====================================================================="

# 1. AUTO-DETECT CLIENT FILES (Idempotency Fix)
REF_FASTA=$(find data/reference . -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | head -n 1)
REF_GTF=$(find data/reference . -maxdepth 1 -type f -name "*.gtf" 2>/dev/null | head -n 1)

if [ -z "$REF_FASTA" ] || [ -z "$REF_GTF" ]; then
    echo "ERROR: Could not find a .fa/.fasta or .gtf file in this folder!"
    exit 1
fi

# 2. SETUP DIRECTORIES & MOVE FILES
mkdir -p data/{raw_fastq,reference} results/{trimmed,aligned,fastqc_reports} R_analysis
mv "$REF_FASTA" data/reference/ 2>/dev/null || true
mv "$REF_GTF" data/reference/ 2>/dev/null || true
mv *.fastq *.fq *.fastq.gz *.fq.gz data/raw_fastq/ 2>/dev/null || true

FASTA_CLEAN=$(basename "$REF_FASTA")
GTF_CLEAN=$(basename "$REF_GTF")
INDEX_BASE="${FASTA_CLEAN%.*}_index"

# 3. BUILD INDEX
if [ ! -f "data/reference/${INDEX_BASE}.1.ht2" ]; then
    echo "==> Building HISAT2 Index (This may take a while)..."
    hisat2-build "data/reference/${FASTA_CLEAN}" "data/reference/${INDEX_BASE}"
fi

# 4. BATCH PROCESS PAIRED-END SAMPLES
shopt -s nullglob
FASTQ_FILES=(data/raw_fastq/*_1.fastq* data/raw_fastq/*_R1.fastq*)
if [ ${#FASTQ_FILES[@]} -eq 0 ]; then
    echo "ERROR: No Paired-End FASTQ files found!"
    exit 1
fi

echo -e "id\ttype" > R_analysis/metadata.tsv

for R1 in "${FASTQ_FILES[@]}"; do
    BASENAME=$(basename "$R1")
    SAMPLE=${BASENAME%_1.fastq*}
    SAMPLE=${SAMPLE%_R1.fastq*}

    # Identify R2 mate (supports uncompressed and gzipped)
    if [ -f "data/raw_fastq/${SAMPLE}_2.fastq" ]; then
        R2="data/raw_fastq/${SAMPLE}_2.fastq"
    elif [ -f "data/raw_fastq/${SAMPLE}_2.fastq.gz" ]; then
        R2="data/raw_fastq/${SAMPLE}_2.fastq.gz"
    elif [ -f "data/raw_fastq/${SAMPLE}_R2.fastq" ]; then
        R2="data/raw_fastq/${SAMPLE}_R2.fastq"
    elif [ -f "data/raw_fastq/${SAMPLE}_R2.fastq.gz" ]; then
        R2="data/raw_fastq/${SAMPLE}_R2.fastq.gz"
    else
        echo "ERROR: Could not find Read 2 for $SAMPLE!"
        continue
    fi

    echo "==> PROCESSING SAMPLE: $SAMPLE"
    echo -e "${SAMPLE}\tnormal" >> R_analysis/metadata.tsv

    fastp -i "$R1" -I "$R2" \
          -o "results/trimmed/${SAMPLE}_1_trimmed.fastq.gz" -O "results/trimmed/${SAMPLE}_2_trimmed.fastq.gz" \
          -h "results/fastqc_reports/fastp_${SAMPLE}_report.html" 2> /dev/null

    hisat2 -x "data/reference/${INDEX_BASE}" \
           -1 "results/trimmed/${SAMPLE}_1_trimmed.fastq.gz" -2 "results/trimmed/${SAMPLE}_2_trimmed.fastq.gz" \
           -S "results/aligned/${SAMPLE}.sam" 2> /dev/null

    samtools sort -@ 4 -o "results/aligned/${SAMPLE}.sorted.bam" "results/aligned/${SAMPLE}.sam"
    rm "results/aligned/${SAMPLE}.sam"
done

# 5. GENERATE COUNT MATRIX
echo "==> Compiling master count matrix..."
featureCounts -T 4 -p -a "data/reference/${GTF_CLEAN}" \
              -o R_analysis/counts_raw.txt results/aligned/*.sorted.bam

echo "==> Formatting matrix for R..."
tail -n +2 R_analysis/counts_raw.txt | cut -f1,7- > R_analysis/counts.tsv
sed -i 's|results/aligned/||g' R_analysis/counts.tsv
sed -i 's|\.sorted\.bam||g' R_analysis/counts.tsv

echo "==========================================================================="
echo "BASH PIPELINE COMPLETE! Edit R_analysis/metadata.tsv before running R."
echo "==========================================================================="
