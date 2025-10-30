### cnvpytor find CNV and anno
source activate cnvpytor
cnvpytor -root file.chr17.pytor -rd /cluster/home/xiejunhao/analysis/20251014-byshiyan/47201/_sam/47201.bqsr.hg38.bam -chrom chr17 -T /genetics/home/stu_liujiyuan/res//hg38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
cnvpytor -root file.chr17.pytor -his 1000 10000 100000
cnvpytor -root file.chr17.pytor -partition 1000 10000 100000
#test_path:/genetics/home/shiyan/tmp/CNV2
# 1) 1 kb 分辨率：更适合检测小事件（≥1–2 kb），但噪声更大
cnvpytor -root file.chr17.pytor -view 1000 <<'END'
set Q0_range 0 0.5                          # 过滤低比对质量读段占比(Q0) >50%的事件
set p_range 0 1e-5                          # e-val1 ≤ 1e-5（更严格可到1e-6）
set size_range 1000 inf                     # 过滤 <1 kb 的事件
set print_filename analysis.chr17.1000.filtered.vcf   # 也可给 .vcf 或 .xlsx
print calls
END

# 2) 10 kb 分辨率：全能型，适合多数体细胞/生殖系的中大事件
cnvpytor -root file.chr17.pytor -view 10000 <<'END'
set Q0_range 0 0.5
set p_range 0 1e-5
set size_range 5000 inf                     # 10 kb 下建议至少 5 kb
set print_filename analysis.chr17.10000.filtered.vcf
print calls
END

# 3) 100 kb 分辨率：稳定捕获 Mb 级大片段，抑制噪声
cnvpytor -root file.chr17.pytor -view 100000 <<'END'
set Q0_range 0 0.5
set p_range 0 1e-5
set size_range 50000 inf                    # 大事件过滤阈值更高
set print_filename analysis.chr17.100000.filtered.vcf
print calls
END

bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\n' analysis.chr17.1000.filtered.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$4,"binsize=10k"}' > calls.1k.bed
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\n' analysis.chr17.10000.filtered.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$4,"binsize=10k"}' > calls.10k.bed
bcftools query -f '%CHROM\t%POS\t%INFO/END\t%SVTYPE\n' analysis.chr17.100000.filtered.vcf | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$3,$4,"binsize=10k"}' > calls.100k.bed

bedtools intersect -a calls.1k.bed  -b calls.10k.bed  -f 0.5 -r -wa -wb | awk '$4==$9' > overlap.1k_10k.bed
bedtools intersect -a calls.10k.bed -b calls.100k.bed -f 0.5 -r -wa -wb | awk '$4==$9' > overlap.10k_100k.bed

# 把出现于任意两组交集里的事件坐标合并去重，作为“共识集”
cat overlap.1k_10k.bed overlap.10k_100k.bed | \
awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4}' | sort -k1,1 -k2,2n | bedtools merge -i - -c 4 -o distinct > consensus.cnv.bed

bedtools intersect -a consensus.cnv.bed -b hg38.genes.bed -wa -wb \
| awk 'BEGIN{OFS="\t"}{
  cS=$2;cE=$3;gS=$6;gE=$7;rel="overlap";
  if (gS>=cS && gE<=cE) rel="inside";
  else if (gS<=cS && gE>=cE) rel="cover";
  else if (gS<cS && gE>=cS && gE<cE) rel="left_bp";
  else if (gS>cS && gS<=cE && gE>cE) rel="right_bp";
  print $1,$2,$3,$4,$5,$6,$7,$8,rel
}' > consensus.cnv.gene.overlap.tsv

# 若没有任何重叠，给出最近基因（及距离）
bedtools closest -a consensus.cnv.bed -b hg38.genes.bed -d \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,"nearest",$9}' \
> consensus.cnv.gene.nearest.tsv
