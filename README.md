# NGS_pipeline

### Bulk RNA-seq

`readlink -f /path/*/*_1.clean.fq.gz|sed 's/_1.clean.fq.gz//g'|awk -F "/" '{print $NF"\t"$0"_1.clean.fq.gz""\t"$0"_2.clean.fq.gz"}' > samples.tsv`

`bash bulk_rna_pipeline.sh samples.tsv SDD`
