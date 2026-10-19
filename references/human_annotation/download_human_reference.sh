#!/usr/bin/env bash
# ==============================================================================
# Script: Download Human Reference Genome (GRCh38) and Gene Annotations (GTF)
# Source: Ensembl Release 110 / GENCODE v44
# Author: Abdul Rehman
# ==============================================================================

set -euo pipefail

DEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${DEST_DIR}"

echo "=== Human Reference Genome & Annotation Downloader ==="
echo "Target Directory: ${DEST_DIR}"

# 1. Official Ensembl GRCh38 Primary Assembly Genome FASTA (~840 MB compressed, ~3.1 GB uncompressed)
GENOME_URL="https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

# 2. Official Ensembl GRCh38 Gene Annotation GTF (~55 MB compressed, ~1.5 GB uncompressed)
GTF_URL="https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

echo "Downloading Human Genomic Annotation GTF..."
curl -L -C - -o "Homo_sapiens.GRCh38.110.gtf.gz" "${GTF_URL}"
echo "Decompressing GTF..."
gzip -dkf "Homo_sapiens.GRCh38.110.gtf.gz"

echo "Downloading Human Primary Assembly FASTA..."
curl -L -C - -o "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" "${GENOME_URL}"
echo "Decompressing FASTA..."
gzip -dkf "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

echo "Reference download complete!"
