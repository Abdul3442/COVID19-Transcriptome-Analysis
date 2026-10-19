#!/usr/bin/env bash
# ==============================================================================
# Step 1: Raw Read Quality Control (FastQC) & Quality Trimming (fastp)
# Author: Abdul Rehman
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_FASTQ_DIR="${PROJECT_ROOT}/data/raw_fastq"
TRIMMED_FASTQ_DIR="${PROJECT_ROOT}/data/trimmed_fastq"
QC_DIR="${PROJECT_ROOT}/results/01_qc"
FASTQC_RAW_DIR="${QC_DIR}/fastqc_raw"
FASTQC_CLEAN_DIR="${QC_DIR}/fastqc_clean"
FASTP_DIR="${QC_DIR}/fastp"

CONDA_BIN="/home/abdul/miniconda3/envs/rnaseq_pipeline/bin"
SRA_BIN="/home/abdul/miniconda3/bin"
export PATH="${CONDA_BIN}:${SRA_BIN}:$PATH"

mkdir -p "${RAW_FASTQ_DIR}" "${TRIMMED_FASTQ_DIR}" "${FASTQC_RAW_DIR}" "${FASTQC_CLEAN_DIR}" "${FASTP_DIR}"

echo -e "\n=== [STEP 1] Raw Read Acquisition & Quality Control ==="

# Representative Biological Replicates from BioProject PRJNA615032 (GSE147507)
# SRR11573904: SARS-CoV-2 Infected Rep 1
# SRR11573908: SARS-CoV-2 Infected Rep 2
# SRR11573892: Mock Control Rep 1
# SRR11573896: Mock Control Rep 2
SAMPLES=("SRR11573904" "SRR11573908" "SRR11573892" "SRR11573896")

for SRA_ID in "${SAMPLES[@]}"; do
    RAW_FASTQ="${RAW_FASTQ_DIR}/${SRA_ID}_1.fastq.gz"
    if [ ! -f "${RAW_FASTQ}" ]; then
        echo ">>> Downloading raw sequencing reads from NCBI SRA for ${SRA_ID}..."
        fastq-dump -X 25000 --split-files --gzip "${SRA_ID}" -O "${RAW_FASTQ_DIR}"
    else
        echo ">>> Raw reads already present for ${SRA_ID}: ${RAW_FASTQ}"
    fi

    # 1. FastQC on Raw Reads
    echo ">>> Running FastQC on raw reads for ${SRA_ID}..."
    fastqc --outdir "${FASTQC_RAW_DIR}" -t 2 "${RAW_FASTQ}"

    # 2. Quality Trimming & Adapter Removal via fastp
    CLEAN_FASTQ="${TRIMMED_FASTQ_DIR}/${SRA_ID}_clean.fastq.gz"
    echo ">>> Running fastp quality & adapter trimming for ${SRA_ID}..."
    fastp \
        -i "${RAW_FASTQ}" \
        -o "${CLEAN_FASTQ}" \
        -h "${FASTP_DIR}/${SRA_ID}_fastp.html" \
        -j "${FASTP_DIR}/${SRA_ID}_fastp.json" \
        -q 20 \
        -u 30 \
        -l 35 \
        --detect_adapter_for_pe \
        -w 2

    # 3. FastQC on Cleaned Reads
    echo ">>> Running post-trimming FastQC on cleaned reads for ${SRA_ID}..."
    fastqc --outdir "${FASTQC_CLEAN_DIR}" -t 2 "${CLEAN_FASTQ}"
done

echo -e "\n>>> Step 1 (Quality Control & Trimming) completed successfully!"
