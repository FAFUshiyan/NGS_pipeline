# cutadapt-crispresso2-batch-pipeline

A lightweight, reproducible bash pipeline for amplicon-seq demultiplexing and CRISPResso2 batch analysis.

This repository provides a one-command workflow to:
1) Demultiplex paired-end FASTQ files by 5' sample barcodes/primers using Cutadapt (named adapters + `{name}` output template),
2) Generate a CRISPResso2 batch settings table,
3) Run CRISPRessoBatch to quantify editing outcomes per sample and produce batch-level comparison summaries when amplicon/guide are shared across samples.

## Input
A tab-separated targets file (e.g., `in-vitro.A9.tsv`) with 4 columns:
1. `sample` (e.g., YKpr1009)
2. `barcode_or_primer` (for Cutadapt demultiplexing)
3. `guide_seq`
4. `amplicon_seq`

## Outputs
- `demux/`: demultiplexed FASTQ files per sample
- `CRISPResso_out_*`: CRISPResso2 outputs per sample plus batch comparison plots/tables

## Dependencies
- Cutadapt
- CRISPResso2 (CRISPRessoBatch)
- GNU coreutils + awk (bash environment)

## Example
./run_demux_crispresso_batch.sh in-vivo.aF9R5-RGA.tsv in-vivo.aF9R5-RGA_OUT M-Edit-LHM24273_L1_1.fq.gz M-Edit-LHM24273_L1_2.fq.gz 12 8
