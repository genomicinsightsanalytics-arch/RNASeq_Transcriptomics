#!/bin/bash
set -e

echo "====================================================================="
echo " UNIVERSAL SINGLE-END RNA-SEQ PIPELINE"
echo "====================================================================="

REF_FASTA=$(find data/reference . -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) 2>/dev/null | head -n 1)
REF_GTF=$(find data/reference . -maxdepth 1 -type f -name "*.gtf" 2>/dev/null | head -n 1)

if [ -z "$REF_FASTA" ] || [ -z "$REF_GTF" ]; then
    echo "ERROR: Could not find a .fa/.fasta or .gtf file in this folder!"
    exit 1
fi

mkdir -p data/{raw_fastq,reference} results/{trimmed,aligned,fastqc_reports} R_analysis
mv "$REF_FASTA" data/reference/ 2>/dev/null || true
mv "$REF_GTF" data/reference/ 2>/dev/null || true
mv *.fastq *.fq *.fastq.gz *.fq.gz data/raw_fastq/ 2>/dev/null || true

FASTA_CLEAN=$(basename "$REF_FASTA")
GTF_CLEAN=$(basename "$REF_GTF")
INDEX_BASE="${FASTA_CLEAN%.*}_index"

if [ ! -f "data/reference/${INDEX_BASE}.1.ht2" ]; then
    echo "==> Building HISAT2 Index..."
    hisat2-build "data/reference/${FASTA_CLEAN}" "data/reference/${INDEX_BASE}"
fi

echo -e "id\ttype" > R_analysis/metadata.tsv

shopt -s nullglob
FASTQ_FILES=(data/raw_fastq/*.fastq*)
for READ in "${FASTQ_FILES[@]}"; do
    BASENAME=$(basename "$READ")
    SAMPLE=${BASENAME%.fastq*}
    SAMPLE=${SAMPLE%.fq*}

    echo "==> PROCESSING SAMPLE: $SAMPLE"
    echo -e "${SAMPLE}\tnormal" >> R_analysis/metadata.tsv

    fastp -i "$READ" -o "results/trimmed/${SAMPLE}_trimmed.fastq.gz" \
          -h "results/fastqc_reports/fastp_${SAMPLE}_report.html" 2> /dev/null

    hisat2 -x "data/reference/${INDEX_BASE}" -U "results/trimmed/${SAMPLE}_trimmed.fastq.gz" \
           -S "results/aligned/${SAMPLE}.sam" 2> /dev/null

    samtools sort -@ 4 -o "results/aligned/${SAMPLE}.sorted.bam" "results/aligned/${SAMPLE}.sam"
    rm "results/aligned/${SAMPLE}.sam"
done

echo "==> Compiling master count matrix..."
featureCounts -T 4 -a "data/reference/${GTF_CLEAN}" -o R_analysis/counts_raw.txt results/aligned/*.sorted.bam

tail -n +2 R_analysis/counts_raw.txt | cut -f1,7- > R_analysis/counts.tsv
sed -i 's|results/aligned/||g' R_analysis/counts.tsv
sed -i 's|\.sorted\.bam||g' R_analysis/counts.tsv

echo "==========================================================================="
echo "BASH PIPELINE COMPLETE! Edit R_analysis/metadata.tsv before running R."
echo "==========================================================================="
