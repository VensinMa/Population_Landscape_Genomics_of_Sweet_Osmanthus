# 步骤 1：提取并格式化基因 ID
awk '{print substr($2, 6)}' LYG.hic.snp.annovar.variant_function > genes_1.txt
awk '{print substr($1, 5, length($1)-5)}' LYG.genomic.pep.Nr5.uniq.annotations > genes_2.txt

# 步骤 2：排序基因 ID 文件
sort genes_1.txt -o genes_1_sorted.txt
sort genes_2.txt -o genes_2_sorted.txt

# 步骤 3：使用 join 命令匹配基因 ID
join -1 1 -2 1 genes_1_sorted.txt genes_2_sorted.txt > matched_genes.txt

# 步骤 4：提取匹配的基因信息
awk 'NR==FNR {a[$1]; next} $2 in a' matched_genes.txt LYG.hic.snp.annovar.variant_function > final_merged_1.txt
awk 'NR==FNR {a[$1]; next} $1 in a' matched_genes.txt LYG.genomic.pep.Nr5.uniq.annotations > final_merged_2.txt

# 步骤 5：合并两个文件
paste final_merged_1.txt final_merged_2.txt > final_merged_data.txt
