# -*- coding: utf-8 -*-

import csv

# SNP ID 文件和目标文件路径
snpid_file = "/home/vensin/workspace/Annovar/result/1766snps.id"
annovar_file = "/home/vensin/workspace/Annovar/result/LYG.hic.snp.nr_annovar.ann"
output_file = "/home/vensin/workspace/Annovar/result/1766_Adaptive_SNP_nr_annovar.ann"

# 读取 SNP IDs 到一个集合中，并且转换为类似 'Superscaffold1:514882' 格式
with open(snpid_file, "r") as snp_file:
    snp_ids = set(line.strip() for line in snp_file)

# 打开目标文件并进行过滤
with open(annovar_file, "r") as ann_file, open(output_file, "w") as out_file:
    reader = csv.reader(ann_file, delimiter='\t')
    writer = csv.writer(out_file, delimiter='\t')
    
    # 遍历每行数据，检查组合的 SNP ID 是否在列表中
    for row in reader:
        scaffold = row[2]  # 假设 SNP 所在的列是第3列
        position = row[3]  # 假设位置是第4列
        snp_id_combined = f"{scaffold}:{position}"  # 拼接为 'Superscaffold1:514882'
        
        # 如果拼接后的 SNP ID 在 SNP ID 列表中，写入该行
        if snp_id_combined in snp_ids:
            writer.writerow(row)

print(f"Extracted lines saved to {output_file}")
