~/bin/Anno.vcf2tsv.top.sh CMTout.annot.vcf.norm.gz CMTout.annot2.tsv

python rare_filter_clinical.py --input CMTout.annot2.tsv --out-prefix CMT.AD --mode AD --maf 0.001 --protein-coding-only

python rare_filter_clinical.py --input CMTout.annot2.tsv --out-prefix CMT.AR --mode AR --maf 0.01 --protein-coding-only --emit-comphet
