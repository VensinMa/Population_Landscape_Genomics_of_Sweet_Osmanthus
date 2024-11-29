# 步骤 1：提取并格式化基因 ID
awk '{print substr($2, 6)}' LYG.hic.snp.annovar.variant_function > result/genes_1.txt
awk '{print substr($1, 5, length($1)-5)}' LYG.genomic.pep.Nr5.uniq.annotations > result/genes_2.txt

# 步骤 2：排序基因 ID 文件
sort result/genes_1.txt -o result/genes_1_sorted.txt
sort result/genes_2.txt -o result/genes_2_sorted.txt

# 步骤 3：使用 join 命令匹配基因 ID
join -1 1 -2 1 result/genes_1_sorted.txt result/genes_2_sorted.txt > result/matched_genes.txt

# 步骤 4：提取匹配的基因信息
awk 'NR==FNR {a[$1]; next} $2 in a' result/matched_genes.txt LYG.hic.snp.annovar.variant_function > result/final_merged_1.txt
awk 'NR==FNR {a[$1]; next} $1 in a' result/matched_genes.txt LYG.genomic.pep.Nr5.uniq.annotations > result/final_merged_2.txt

# 步骤 5：合并两个文件
paste result/final_merged_1.txt result/final_merged_2.txt > result/final_merged_data.txt
