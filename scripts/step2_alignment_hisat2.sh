#!/usr/bin/env bash
# ==============================================================================
# Step 2: Read Alignment & Mapping (HISAT2) and BAM Processing (samtools)
# Author: Abdul Rehman
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TRIMMED_FASTQ_DIR="${PROJECT_ROOT}/data/trimmed_fastq"
REF_DIR="${PROJECT_ROOT}/references/viral_genome"
ALIGN_DIR="${PROJECT_ROOT}/results/02_alignment"
INDEX_DIR="${PROJECT_ROOT}/references/hisat2_index"

CONDA_BIN="/home/abdul/miniconda3/envs/rnaseq_pipeline/bin"
export PATH="${CONDA_BIN}:$PATH"

mkdir -p "${ALIGN_DIR}" "${INDEX_DIR}"

REF_FASTA="${REF_DIR}/SARS_CoV_2_Wuhan_Hu_1.fasta"
INDEX_PREFIX="${INDEX_DIR}/sars_cov_2_index"

echo -e "\n=== [STEP 2] Read Alignment with HISAT2 & SAMtools ==="

# 1. Build HISAT2 Genome Index if not already present
if [ ! -f "${INDEX_PREFIX}.1.ht2" ]; then
    echo ">>> Building HISAT2 index for ${REF_FASTA}..."
    hisat2-build -p 2 "${REF_FASTA}" "${INDEX_PREFIX}"
else
    echo ">>> HISAT2 index already built at ${INDEX_PREFIX}"
fi

# 2. Align Clean Reads and Process with samtools
SAMPLES=("SRR11573904" "SRR11573908" "SRR11573892" "SRR11573896")

for SRA_ID in "${SAMPLES[@]}"; do
    CLEAN_FASTQ="${TRIMMED_FASTQ_DIR}/${SRA_ID}_clean.fastq.gz"
    SAM_FILE="${ALIGN_DIR}/${SRA_ID}.sam"
    SORTED_BAM="${ALIGN_DIR}/${SRA_ID}_sorted.bam"
    FLAGSTAT_FILE="${ALIGN_DIR}/${SRA_ID}_flagstat.txt"
    SUMMARY_LOG="${ALIGN_DIR}/${SRA_ID}_hisat2_summary.txt"

    echo ">>> Aligning ${SRA_ID} with HISAT2..."
    hisat2 \
        -p 2 \
        -x "${INDEX_PREFIX}" \
        -U "${CLEAN_FASTQ}" \
        -S "${SAM_FILE}" \
        --summary-file "${SUMMARY_LOG}"

    echo ">>> Sorting & indexing BAM with samtools for ${SRA_ID}..."
    samtools view -bS -@ 2 "${SAM_FILE}" | samtools sort -@ 2 -o "${SORTED_BAM}" -
    samtools index "${SORTED_BAM}"

    echo ">>> Generating alignment flagstat report for ${SRA_ID}..."
    samtools flagstat "${SORTED_BAM}" > "${FLAGSTAT_FILE}"

    # Clean intermediate SAM
    rm -f "${SAM_FILE}"
    echo ">>> Alignment completed for ${SRA_ID}: ${SORTED_BAM}"
done

echo -e "\n>>> Step 2 (HISAT2 Alignment & BAM Processing) completed successfully!"
