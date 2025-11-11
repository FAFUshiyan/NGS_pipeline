###计算每个区域的平均覆盖深度
mosdepth --by targets.bed --threads 8 S01 S01.bam
###保留平均深度≥10×的目标区段
for k in `cat  controll.txt`; do zcat $k.regions.bed.gz | awk '$4>=10{print $1"\t"$2"\t"$3}' > $k.mean_ge10x.bed; done
# 多样本：统计每个完全相同的区段，被≥多少个样本满足“平均≥10×”, 62个样本的90%≈56
cat *.mean_ge10x.bed \
| sort -k1,1 -k2,2n -k3,3n \
| awk 'BEGIN{OFS="\t"}{print $1,$2,$3,1}' \
| bedtools groupby -g 1,2,3 -c 4 -o sum \
| awk -v n=56 '$4>=n{print $1"\t"$2"\t"$3}' \
> callable.mean10x.90p.bed
### 保留“共同可调用区”的变异
bcftools view -R ../callable.mean10x.90p.bed -Oz -o merged.65samples.masked.vcf.gz merged.65samples.vcf.gz
### vcf norm
bcftools norm --threads 8 -f ~/ref_data/Human/DNA/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna -m -both -Oz -o merged.65samples.masked.norm.vcf.gz merged.65samples.masked.vcf.gz
### 保留缺失率小于20%的位点
~/miniforge3/envs/vcftools/bin/vcftools --gzvcf merged.65samples.masked.norm.vcf.gz   --max-missing 0.8   --recode --recode-INFO-all  --out filtered
### 保留注释top的内容
cat filtered.recode.vcf |~/bin/snpEff/scripts/_OLD/vcfAnnFirst.py  > merged.65samples.norm.top.vcf
### 保留SNP
bcftools view -v snps -Oz -o merged.65samples.norm.snp.top.vcf.gz merged.65samples.norm.top.vcf
### 提取注释中的关键内容，为下一步分析提供基础
zsh extract.sh merged.65samples.norm.snp.top.vcf.gz merged.65samples.norm.snp.top.tsv
### 生成pheno_covars.tsv 10代表PCA的数目
zsh make_pheno_covars.sh merged.65samples.norm.snp.top.vcf.gz output.65sample 10 case.ID control.ID
### 将变异分为三个常用面具（LoF / Damaging missense / Splice）
python make_setid.py merged.65samples.norm.snp.top.tsv  --bim output.65sample.bed1u.bim  --maf 0.0001 --revel 0.5 --cadd 20 --alphamiss 0.5 --spliceai 0.5 --id-style colon --force-chr no
### 执行SKAT分析
/public/software/VersionHub/R/4.4.2/Rscript SKAT.R
