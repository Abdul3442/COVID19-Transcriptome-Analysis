#!/usr/bin/env bash
# ==============================================================================
# Step 3: Read Quantification (featureCounts) & Expression Matrix (TPM/FPKM)
# Author: Abdul Rehman
# ==============================================================================

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ALIGN_DIR="${PROJECT_ROOT}/results/02_alignment"
REF_DIR="${PROJECT_ROOT}/references/viral_genome"
QUANT_DIR="${PROJECT_ROOT}/results/03_quantification"

CONDA_BIN="/home/abdul/miniconda3/envs/rnaseq_pipeline/bin"
PYTHON_BIN="${PROJECT_ROOT}/venv/bin/python3"
export PATH="${CONDA_BIN}:$PATH"

mkdir -p "${QUANT_DIR}"

ANNOTATION_GFF="${REF_DIR}/SARS_CoV_2_ASM985889v3_genomic.gff"
COUNTS_OUT="${QUANT_DIR}/featurecounts_raw.txt"

echo -e "\n=== [STEP 3] Gene Quantification with featureCounts ==="

BAM_FILES=(
    "${ALIGN_DIR}/SRR11573904_sorted.bam"
    "${ALIGN_DIR}/SRR11573908_sorted.bam"
    "${ALIGN_DIR}/SRR11573892_sorted.bam"
    "${ALIGN_DIR}/SRR11573896_sorted.bam"
)

echo ">>> Running featureCounts across sorted BAM files..."
featureCounts \
    -a "${ANNOTATION_GFF}" \
    -t gene \
    -g ID \
    -o "${COUNTS_OUT}" \
    -T 2 \
    "${BAM_FILES[@]}"

echo ">>> Formatting counts table and computing TPM matrix via Python..."
"${PYTHON_BIN}" -c "
import pandas as pd
import numpy as np

# Load featureCounts output (skipping comment line 1)
df = pd.read_csv('${COUNTS_OUT}', sep='\t', skiprows=1)

# Clean column headers
cols = ['Geneid', 'Length']
sample_cols = [c for c in df.columns if c.endswith('_sorted.bam')]
clean_sample_names = [c.split('/')[-1].replace('_sorted.bam', '') for c in sample_cols]

count_df = df[['Geneid', 'Length'] + sample_cols].copy()
count_df.columns = ['gene_id', 'gene_length_bp'] + clean_sample_names

# Export Clean Raw Counts
raw_counts_matrix = count_df.set_index('gene_id')
raw_matrix_file = '${QUANT_DIR}/featurecounts_raw_counts.tsv'
raw_counts_matrix[clean_sample_names].to_csv(raw_matrix_file, sep='\t')
print(f'Raw count matrix saved to: {raw_matrix_file}')

# Compute Transcripts Per Million (TPM)
length_kb = count_df['gene_length_bp'].values / 1000.0
rpk = count_df[clean_sample_names].div(length_kb, axis=0)
scaling_factor = rpk.sum(axis=0) / 1e6
scaling_factor[scaling_factor == 0] = 1.0
tpm = rpk.div(scaling_factor, axis=1)
tpm.index = count_df['gene_id']

tpm_file = '${QUANT_DIR}/featurecounts_tpm_matrix.tsv'
tpm.to_csv(tpm_file, sep='\t')
print(f'TPM expression matrix saved to: {tpm_file}')

# Copy summary
summary_file = '${COUNTS_OUT}.summary'
clean_summary = pd.read_csv(summary_file, sep='\t')
clean_summary.to_csv('${QUANT_DIR}/featurecounts_summary.tsv', sep='\t', index=False)
"

echo -e "\n>>> Step 3 (Quantification & TPM Matrix) completed successfully!"
