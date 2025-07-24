# NGS_pipeline

## Overview

This repository provides a one‑stop **Bulk RNA‑seq** processing pipeline—from raw FASTQ through alignment, gene counting, and TPM normalization—implemented in a single shell script (`bulk_rna_pipeline.sh`). You supply:  
1. A **samples.tsv** file listing your sample IDs and FASTQ paths  
2. An **output prefix** for counts/TPM files  

The script then:  
- Runs **STAR** to align and produce sorted BAMs  
- Uses **featureCounts** to generate a raw gene‑level count matrix  
- Applies **rnanorm tpm** to compute a TPM matrix  

## Contents
├── bulk_rna_pipeline.sh # All‑in‑one pipeline driver

├── samples.tsv # Sample sheet (Tab‑separated)

└── README.md # This file
## samples.tsv

Each line is a sample, with **three** tab‑separated columns:

Sample1 /path/to/Sample1_1.clean.fq.gz /path/to/Sample1_2.clean.fq.gz

Sample2 /path/to/Sample2_1.clean.fq.gz /path/to/Sample2_2.clean.fq.gz

## Usage

```bash
bash bulk_rna_pipeline.sh samples.tsv <output_prefix>
samples.tsv    — your sample sheet as above

<output_prefix> — prefix for the counts & TPM files (e.g. MyProject)

## Example

bash bulk_rna_pipeline.sh samples.tsv SDD