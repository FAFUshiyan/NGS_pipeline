out="CRISPResso_editing_summary.tsv"
echo -e "Sample\tModified%\tReads_in_input\tReads_aligned_all_amplicons" > "$out"

for f in ./CRISPResso_on_YKpr*/CRISPResso_quantification_of_editing_frequency.txt; do
  sample=$(basename "$(dirname "$f")" | sed 's/^CRISPResso_on_//')
  awk -F'\t' -v s="$sample" '
    NR==1 {next}                         # 跳过表头
    $0 ~ /^[[:space:]]*$/ {next}         # 跳过空行
    { gsub(/\r$/, "", $0) }              # 兼容CRLF
    { print s "\t" $3 "\t" $4 "\t" $5; exit }   # 取第一条数据行(通常就是Reference那行)
  ' "$f" >> "$out"
done
