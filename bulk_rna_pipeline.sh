#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# bulk_rna_pipeline.sh
#
# 一体化从 FASTQ 到排序 BAM、基因计数、TPM 矩阵输出的完整流水线
#
# 用法：
#   ./bulk_rna_pipeline.sh <samples.tsv> <output_prefix>
#
#   samples.tsv     三列 Tab：SampleID    R1.fastq.gz    R2.fastq.gz
#   output_prefix   用于生成 counts 和 TPM 文件的前缀
# ==============================================================================

### —— 配置区 —— ###
GENOME_DIR="${HOME}/data/RNA-seq/GRCh38_STAR"
GTF_FILE="${HOME}/data/RNA-seq/Homo_sapiens.GRCh38.108.gtf"
THREADS=12

### —— 参数校验 —— ###
if [[ $# -ne 2 ]]; then
  echo "Usage: $0 <samples.tsv> <output_prefix>" >&2
  exit 1
fi
SAMPLES_TSV="$1"
PREFIX="$2"

if [[ ! -f "${SAMPLES_TSV}" ]]; then
  echo "Error: 找不到样本表 ${SAMPLES_TSV}" >&2
  exit 1
fi

echo "=== Bulk RNA‑seq Pipeline ==="
echo "Samples:   ${SAMPLES_TSV}"
echo "Prefix:    ${PREFIX}"
echo "Threads:   ${THREADS}"
echo

# ------------------------------------------------------------------------------
# 1. STAR 比对并生成排序 BAM
# ------------------------------------------------------------------------------
echo ">> Step 1/3: STAR alignment"
while IFS=$'\t' read -r SAMPLE R1 R2; do
  echo "--> ${SAMPLE}"
  STAR \
    --genomeDir "${GENOME_DIR}" \
    --runThreadN "${THREADS}" \
    --readFilesIn "${R1}" "${R2}" \
    --readFilesCommand zcat \
    --outFileNamePrefix "${SAMPLE}." \
    --outSAMtype BAM SortedByCoordinate \
    --outBAMsortingThreadN "${THREADS}" \
    --sjdbGTFfile "${GTF_FILE}"
done < "${SAMPLES_TSV}"
echo

# ------------------------------------------------------------------------------
# 2. featureCounts 统计基因计数
# ------------------------------------------------------------------------------
echo ">> Step 2/3: featureCounts"
featureCounts \
  -T "${THREADS}" \
  -t exon \
  -g gene_id \
  -Q 10 \
  --primary \
  -p \
  -a "${GTF_FILE}" \
  -F GTF \
  -o "${PREFIX}.counts" \
  ./*.bam
echo "Generated counts: ${PREFIX}.counts.txt"
echo

# ------------------------------------------------------------------------------
# 3. 生成 TPM 矩阵
# ------------------------------------------------------------------------------
echo ">> Step 3/3: generate TPM matrix"
bash ~/bin/get.tpm.sh ${PREFIX}.counts
echo "Generated TPM matrix: ${PREFIX}.tpm.tsv"
echo

echo "=== Pipeline done ==="
