java -jar ~/bin/snpEff/SnpSift.jar extractFields -s "," -e "." "$1" \
  CHROM POS REF ALT ID \
  "ANN[0].GENE" "ANN[0].FEATUREID" "ANN[0].EFFECT" "ANN[0].IMPACT" \
  "ANN[0].BIOTYPE" "ANN[0].RANK" "ANN[0].HGVS_C" "ANN[0].HGVS_P" \
  gnomAD_cadd_phred REVEL_max SpliceAI_DS_max am_pathogenicity am_class \
  gnomAD_AF_afr gnomAD_AF_amr gnomAD_AF_asj gnomAD_AF_eas \
  gnomAD_AF_fin gnomAD_AF_mid gnomAD_AF_nfe gnomAD_AF_sas \
  WBBC_AF WBBC_AF_Central WBBC_AF_Lingnan WBBC_AF_North WBBC_AF_South \
  "GEN[0].GT" "GEN[0].DP" "GEN[0].AD" "GEN[0].GQ" \
  > "$2"
