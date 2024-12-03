# -*- coding: utf-8 -*-

import csv

# 读取nr注释和annovar注释文件
nr_annotations_file = "/home/vensin/workspace/nr.annotations/LYG.longest.pep.uniq.Nr.modified.annotations"
annovar_annotations_file = "/home/vensin/workspace/Annovar/LYG.hic.snp.annovar.variant_function"
output_file = "/home/vensin/workspace/Annovar/result/merged_output.txt"

# 存储第一个文件数据，按基因ID索引
annotations_data = {}
with open(nr_annotations_file, "r") as af:
    reader = csv.reader(af, delimiter='\t')
    for row in reader:
        gene_id = row[0]  # 基因ID在第一列
        annotations_data[gene_id] = row

# 读取第二个文件并提取基因ID，合并相关信息
with open(annovar_annotations_file, "r") as vf, open(output_file, "w") as out:
    reader = csv.reader(vf, delimiter='\t')
    writer = csv.writer(out, delimiter='\t')
    
    for row in reader:
        # 假设基因ID嵌入在描述字段中，如 gene-LYG000001，位于第二列
        description = row[1]  
        # 提取基因ID
        if "gene-" in description:
            gene_id = description.split("gene-")[1].split("(")[0]  # 提取基因ID
        else:
            continue
        
        # 如果基因ID在第一个文件中，进行合并
        if gene_id in annotations_data:
            merged_row = annotations_data[gene_id] + row
            writer.writerow(merged_row)

print(f"合并完成，结果保存到 {output_file}")
