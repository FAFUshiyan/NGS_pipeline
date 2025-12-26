#!/usr/bin/env bash
set -euo pipefail

# 用法:
#   bash run_A9_vitro.sh in-vitro.A9.tsv OUTDIR Edit-LHM20712_L1_1.fq.gz Edit-LHM20712_L1_2.fq.gz 12 8
#
# 参数:
#   $1 targets.tsv
#   $2 输出目录
#   $3 R1 fastq.gz
#   $4 R2 fastq.gz
#   $5 cutadapt线程(默认12)
#   $6 CRISPResso线程(默认8)

targets_tsv="$1"
outdir="$2"
r1="$3"
r2="$4"
cutadapt_threads="${5:-12}"
crisp_threads="${6:-8}"

mkdir -p "$outdir"
demux_dir="$outdir/demux"
crisp_dir="$outdir/CRISPResso_out"
mkdir -p "$demux_dir" "$crisp_dir"

# 1) 生成一个 barcodes.fasta（不要用 .tsv 后缀，避免格式识别问题）
#    Cutadapt 支持从 fasta 读取 named adapters，并用 {name} 做 demultiplex。:contentReference[oaicite:2]{index=2}
barcodes_fa="$outdir/barcodes.A9.fasta"
awk -F'\t' '
  NF>=2 && $1!="" && $2!="" {
    gsub(/\r$/, "", $0);   # 去掉可能的CRLF
    print ">" $1 "\n" $2
  }
' "$targets_tsv" > "$barcodes_fa"

# 2) Cutadapt 一次性 demultiplex（输出文件名里必须同时包含 {name} in -o 与 -p 才能 paired-end demux）:contentReference[oaicite:3]{index=3}
cutadapt \
  -e 0 --no-indels \
  -g "^file:${barcodes_fa}" \
  --discard-untrimmed \
  -o "${demux_dir}/demux-vitro.{name}.R1.fq.gz" \
  -p "${demux_dir}/demux-vitro.{name}.R2.fq.gz" \
  -j "$cutadapt_threads" \
  "$r1" "$r2"

# 3) 生成 CRISPRessoBatch 的 batch_settings（表头是参数名，不带 --，官方文档明确此格式）:contentReference[oaicite:4]{index=4}
batch_tsv="$outdir/CRISPResso.batch_settings.tsv"
{
  printf "name\tfastq_r1\tfastq_r2\tamplicon_seq\tguide_seq\n"
  awk -F'\t' -v demux="$demux_dir" '
    NF>=4 && $1!="" {
      gsub(/\r$/, "", $0);
      printf "%s\t%s/demux-vitro.%s.R1.fq.gz\t%s/demux-vitro.%s.R2.fq.gz\t%s\t%s\n", \
             $1, demux, $1, demux, $1, $4, $3
    }
  ' "$targets_tsv"
} > "$batch_tsv"

# 4) 跑 CRISPRessoBatch（-bs 指定 batch_settings，-bo 指定 batch 输出目录）:contentReference[oaicite:5]{index=5}
if command -v CRISPRessoBatch >/dev/null 2>&1; then
  CRISPRessoBatch \
    -bs "$batch_tsv" \
    -bo "$crisp_dir" \
    -p "$crisp_threads" \
    --skip_failed
else
  echo "[ERROR] 未找到 CRISPRessoBatch。请确认你安装的是 CRISPResso2（带 Batch 工具）。" >&2
  exit 1
fi

# 5) 归档 inputs
cp -f "$targets_tsv" "$outdir/"
echo "[DONE] 输出目录: $outdir"
