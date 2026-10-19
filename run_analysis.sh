#!/usr/bin/env bash
# ==============================================================================
# Complete 6-Step End-to-End Bioinformatics RNA-seq Analysis Pipeline
# Language: Bash & Python
# Author: Abdul Rehman
# ==============================================================================

set -euo pipefail

BOLD="\033[1m"
GREEN="\033[0;32m"
BLUE="\033[0;34m"
CYAN="\033[0;36m"
YELLOW="\033[1;33m"
RED="\033[0;31m"
RESET="\033[0m"

PIPELINE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export MPLCONFIGDIR="/tmp/mpl"
mkdir -p "${MPLCONFIGDIR}"

log_header() {
    echo -e "\n${BOLD}${CYAN}======================================================================${RESET}"
    echo -e "${BOLD}${CYAN}   $1${RESET}"
    echo -e "${BOLD}${CYAN}======================================================================${RESET}"
}

log_step() {
    echo -e "\n${BOLD}${BLUE}>>> [STEP $1] $2...${RESET}"
}

log_success() {
    echo -e "${BOLD}${GREEN}[SUCCESS] $1${RESET}"
}

TOTAL_START=$(date +%s)

log_header "Complete 6-Step RNA-seq Transcriptome Analysis Pipeline"
echo -e "${BOLD}Project Directory:${RESET} ${PIPELINE_ROOT}"
echo -e "${BOLD}Investigator:${RESET}      Abdul Rehman"
echo -e "${BOLD}Timestamp:${RESET}         $(date '+%Y-%m-%d %H:%M:%S')"

# ------------------------------------------------------------------------------
# STEP 1: Quality Control (FastQC) & Trimming (fastp)
# ------------------------------------------------------------------------------
log_step "1/6" "Quality Control: FastQC & fastp Adapter/Quality Trimming"
T_START=$(date +%s)
"${PIPELINE_ROOT}/scripts/step1_qc_and_trimming.sh"
T_END=$(date +%s)
log_success "Step 1 completed in $((T_END - T_START))s"

# ------------------------------------------------------------------------------
# STEP 2: Read Alignment / Mapping (HISAT2 & SAMtools)
# ------------------------------------------------------------------------------
log_step "2/6" "Read Alignment: HISAT2 Mapping & SAMtools BAM Sorting/Indexing"
T_START=$(date +%s)
"${PIPELINE_ROOT}/scripts/step2_alignment_hisat2.sh"
T_END=$(date +%s)
log_success "Step 2 completed in $((T_END - T_START))s"

# ------------------------------------------------------------------------------
# STEP 3: Quantification (featureCounts & TPM Calculation)
# ------------------------------------------------------------------------------
log_step "3/6" "Read Quantification: featureCounts & TPM/FPKM Matrix Construction"
T_START=$(date +%s)
"${PIPELINE_ROOT}/scripts/step3_quantification.sh"
T_END=$(date +%s)
log_success "Step 3 completed in $((T_END - T_START))s"

# ------------------------------------------------------------------------------
# STEP 4: Differential Expression Analysis (DESeq2/Empirical Bayes)
# ------------------------------------------------------------------------------
log_step "4/6" "Differential Expression Analysis (DEA): Wald Testing & BH-FDR Correction"
T_START=$(date +%s)
"${PIPELINE_ROOT}/venv/bin/python3" "${PIPELINE_ROOT}/scripts/step4_differential_expression.py"
T_END=$(date +%s)
log_success "Step 4 completed in $((T_END - T_START))s"

# ------------------------------------------------------------------------------
# STEP 5: Functional Annotation & Pathway Enrichment (GO & KEGG)
# ------------------------------------------------------------------------------
log_step "5/6" "Functional Annotation: Gene Ontology (GO) & KEGG Pathway Enrichment"
T_START=$(date +%s)
"${PIPELINE_ROOT}/venv/bin/python3" "${PIPELINE_ROOT}/scripts/step5_pathway_enrichment.py"
T_END=$(date +%s)
log_success "Step 5 completed in $((T_END - T_START))s"

# ------------------------------------------------------------------------------
# STEP 6: Data Visualization & Interpretation (Publication Figures)
# ------------------------------------------------------------------------------
log_step "6/6" "Data Visualization: Publication Figures (PCA, Heatmaps, Volcano, Bubble)"
T_START=$(date +%s)
"${PIPELINE_ROOT}/venv/bin/python3" "${PIPELINE_ROOT}/scripts/step6_generate_figures.py"
T_END=$(date +%s)
log_success "Step 6 completed in $((T_END - T_START))s"

TOTAL_END=$(date +%s)

# ------------------------------------------------------------------------------
# Summary Banner
# ------------------------------------------------------------------------------
log_header "All 6 Pipeline Steps Completed Successfully in $((TOTAL_END - TOTAL_START))s!"
echo -e "${BOLD}${GREEN}Generated Pipeline Assets:${RESET}"
echo -e "  - ${BOLD}Quality Control:${RESET}     ${PIPELINE_ROOT}/results/01_qc/ (Pure text TSV metrics & summary tables)"
echo -e "  - ${BOLD}Alignment & BAMs:${RESET}    ${PIPELINE_ROOT}/results/02_alignment/ (Sorted BAM, BAI, flagstat)"
echo -e "  - ${BOLD}Quantification:${RESET}      ${PIPELINE_ROOT}/results/03_quantification/ (Raw counts & TPM matrix)"
echo -e "  - ${BOLD}Differential Expr:${RESET}   ${PIPELINE_ROOT}/results/04_dea/ (Full statistics & top DEGs)"
echo -e "  - ${BOLD}Pathway Enrichment:${RESET}  ${PIPELINE_ROOT}/results/05_enrichment/ (GO & KEGG tables)"
echo -e "  - ${BOLD}Publication Figures:${RESET} ${PIPELINE_ROOT}/results/06_figures/ (300 DPI vector PNGs)"
echo ""
