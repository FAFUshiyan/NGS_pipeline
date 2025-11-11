#!/usr/bin/env bash
set -euo pipefail

# ===== 参数区（按需修改/传参） =====
VCF="${1:-cohort.exome.vcf.gz}"     # 输入：多样本 VCF/BCF（建议仅外显子区域）
OUT="${2:-cohort}"                  # 中间产物与输出前缀
KPC="${3:-10}"                      # 计算前 K 个主成分
CASE_LIST="${4:-}"                  # 可选：病例样本名单（一列IID）
CTRL_LIST="${5:-}"                  # 可选：对照样本名单（一列IID）

echo "[INFO] VCF=${VCF}"
echo "[INFO] OUT prefix=${OUT}, PCs=${KPC}"

# 0) 样本ID
bcftools query -l "${VCF}" > "${OUT}.samples.id"

# 1) 生成表型（pheno）；若提供 case/control 名单则赋值 1/0，否则先置 NA
#   文件：OUT.pheno.tsv => 列头：#IID pheno
awk 'BEGIN{OFS="\t"; print "#IID","pheno"}{print $1,"NA"}' "${OUT}.samples.id" > "${OUT}.pheno.tsv"

if [[ -n "${CASE_LIST}" && -s "${CASE_LIST}" ]]; then
  awk 'NR==FNR{a[$1]=1;next} BEGIN{OFS="\t"} FNR==1{print $0; next} {p=$2; if($1 in a)p=1; print $1,p}' \
    "${CASE_LIST}" "${OUT}.pheno.tsv" > "${OUT}.pheno.tmp" && mv "${OUT}.pheno.tmp" "${OUT}.pheno.tsv"
fi

if [[ -n "${CTRL_LIST}" && -s "${CTRL_LIST}" ]]; then
  awk 'NR==FNR{a[$1]=1;next} BEGIN{OFS="\t"} FNR==1{print $0; next} {p=$2; if($1 in a)p=0; print $1,p}' \
    "${CTRL_LIST}" "${OUT}.pheno.tsv" > "${OUT}.pheno.tmp" && mv "${OUT}.pheno.tmp" "${OUT}.pheno.tsv"
fi

# 2) VCF → PLINK2 二进制；含X染色体时按人类PAR边界拆分 + 统一变异ID
plink2 --vcf "${VCF}" \
  --psam merge.sex.ID \
  --split-par hg38 \
  --set-all-var-ids @:#:\$r:\$a \
  --new-id-max-allele-len 50 missing \
  --make-bed --out ${OUT}.bed0_pre
# 说明：@=染色体，#=位置，$r/$a=REF/第一ALT；Bash 里 $ 需要转义。超过 400 的等位仅在“ID字符串里”截断，而不是变成“.”。:contentReference[oaicite:3]{index=3}

# 2) 先真正去重，再继续
plink2 --bfile ${OUT}.bed0_pre \
  --rm-dup exclude-all list \
  --make-bed --out ${OUT}.bed0
# 说明：会写出 ${OUT}.bed0.rmdup.list，列出被判定为重复的 ID；这一步才是让 ID 唯一化的关键。:contentReference[oaicite:4]{index=4}

# 4) 保险起见，再检查/去一次重（可选）
plink2 --bfile ${OUT}.bed0 --rm-dup exclude-all --make-bed --out ${OUT}.bed1u

# 5) 现在再做 LD 修剪 + PCA 就不会再报 “requires unique variant IDs”
plink2 --bfile ${OUT}.bed1u --maf 0.05 --indep-pairwise 200 50 0.1 --out ${OUT}.prune
plink2 --bfile ${OUT}.bed1u --extract ${OUT}.prune.prune.in --pca 10 approx --out ${OUT}.pca


# 5) 提取 sex（建议与用于PCA的同一数据集一致；此处用 bed1u）
awk 'BEGIN{OFS="\t"; print "#IID","sex"}{print $2,$5}' "${OUT}.bed1u.fam" > "${OUT}.sex.tsv"

# 6) 合并 pheno + sex （两边都是TAB、以IID为第一列）
{
  printf "%s\n" "IID	pheno	sex"
  tail -n +2 "${OUT}.pheno.tsv" | sort -k1,1
} > "${OUT}.pheno.sorted.tsv"

{
  printf "%s\n" "#IID	sex"
  tail -n +2 "${OUT}.sex.tsv" | sort -k1,1
} > "${OUT}.sex.sorted.tsv"

join -t $'\t' -a 1 -e NA -o auto -1 1 -2 1 \
  <(tail -n +2 "${OUT}.pheno.sorted.tsv") \
  <(tail -n +2 "${OUT}.sex.sorted.tsv") \
  > "${OUT}.ph_se.tsv"

# 7) 将 .eigenvec 转成“TAB分隔、IID首列”的文件，再与 ph_se.tsv 按 IID 合并
#    .eigenvec 是空白分隔，前两列一般为 FID 和 IID；我们取 IID + PCs
awk 'BEGIN{OFS="\t"}
NR==1{
  # 头：去掉前两列（FID IID），保留 PC 列名
  hdr=""; for(i=3;i<=NF;i++) hdr = hdr (i==3? $i:OFS $i);
  print "IID\t" hdr;
  next
}
{
  # 行：第2列是IID，后面接PC值
  printf "%s", $2;
  for(i=3;i<=NF;i++) printf OFS $i;
  printf "\n";
}' "${OUT}.pca.eigenvec" > "${OUT}.evec4join.tsv"

# 准备表头与内容
PC_HEADER=$(head -1 "${OUT}.evec4join.tsv" | cut -f2-)
tail -n +2 "${OUT}.evec4join.tsv" | sort -k1,1 > "${OUT}.eigenvec.sorted"

# 最终合并（左连接，IID为键；未命中的PC用NA）
printf "IID\tpheno\tsex\t%s\n" "${PC_HEADER}" > pheno_covars.tsv
join -t $'\t' -a 1 -e NA -o auto -1 1 -2 1 \
  <(sort -k1,1 "${OUT}.ph_se.tsv") \
  "${OUT}.eigenvec.sorted" \
  >> pheno_covars.tsv

echo "[DONE] pheno_covars.tsv 已生成"

