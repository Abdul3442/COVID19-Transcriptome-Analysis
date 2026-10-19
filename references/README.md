# Reference Genomes and Annotations Guide

This directory contains the reference genomes and gene annotation files used for mapping, transcript quantification, and biological interpretation in this COVID-19 transcriptome analysis project.

---

## 1. Directory Structure

```
references/
├── README.md
├── viral_genome/
│   ├── SARS_CoV_2_Wuhan_Hu_1.fasta         # Complete SARS-CoV-2 reference sequence (NC_045512.2)
│   └── SARS_CoV_2_ASM985889v3_genomic.gff  # Viral genomic annotation (S, N, M, E, ORF1ab, etc.)
└── human_annotation/
    └── download_human_reference.sh         # Shell script to fetch Ensembl GRCh38 FASTA & GTF
```

---

## 2. SARS-CoV-2 Viral Genome

- **Organism**: Severe acute respiratory syndrome coronavirus 2 (SARS-CoV-2)
- **Isolate**: Wuhan-Hu-1
- **NCBI Accession**: [`NC_045512.2`](https://www.ncbi.nlm.nih.gov/nuccore/NC_045512.2) / Assembly `GCF_009858895.2` (`ASM985889v3`)
- **Length**: 29,903 bp (ssRNA)
- **Annotation (`.gff`)**: Contains genomic coordinates for structural proteins (Spike `S`, Nucleocapsid `N`, Membrane `M`, Envelope `E`) and non-structural polyproteins (`ORF1ab`, `ORF3a`, `ORF6`, `ORF7a`, `ORF8`, `ORF10`).
- **File location**: `references/viral_genome/`

---

## 3. Human Host Reference Genome & Annotation (GRCh38)

- **Organism**: *Homo sapiens*
- **Assembly**: Genome Reference Consortium Human Build 38 patch 14 (**GRCh38.p14**)
- **Annotation**: Ensembl Release 110 / GENCODE v32–v44
- **Transcriptome Gene Counts**: The count matrices in `data/raw/` and `data/processed/` are quantified against HGNC gene symbols mapped to the GRCh38 human genome build.
- **Downloading full raw Human FASTA & GTF**: If you wish to re-align raw FASTQ sequencing reads using STAR, HISAT2, or Salmon, execute:
  ```bash
  cd references/human_annotation
  ./download_human_reference.sh
  ```
