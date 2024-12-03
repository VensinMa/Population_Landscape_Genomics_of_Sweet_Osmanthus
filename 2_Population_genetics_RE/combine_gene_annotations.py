# -*- coding: utf-8 -*-

import csv
import time

# 记录脚本开始时间
start_time = time.time()

# 读取nr注释和annovar注释文件
nr_annotations_file = "/home/vensin/workspace/nr.annotations/LYG.longest.pep.uniq.Nr.modified.annotations"
annovar_annotations_file = "/home/vensin/workspace/Annovar/LYG.hic.snp.annovar.variant_function"
output_file = "/home/vensin/workspace/Annovar/result/LYG.hic.snp.nr_annovar.ann"

# 存储第一个文件数据，按基因ID索引
annotations_data = {}
with open(nr_annotations_file, "r") as af:
    reader = csv.reader(af, delimiter='\t')
    for row in reader:
        gene_id = row[0]  # 基因ID在第一列
        annotations_data[gene_id] = row

# 使用生成器来按行处理文件，避免一次性加载整个文件
def process_variant_file(file_path):
    with open(file_path, "r") as vf:
        reader = csv.reader(vf, delimiter='\t')
        for row in reader:
            description = row[1]  # 基因ID在第二列的描述字段
            if "gene-" in description:
                gene_id = description.split("gene-")[1].split("(")[0]  # 提取基因ID
                yield gene_id, row

# 打开输出文件
with open(output_file, "w") as out:
    writer = csv.writer(out, delimiter='\t')
    
    # 处理annovar注释文件，查找匹配的基因ID并合并数据
    for gene_id, row in process_variant_file(annovar_annotations_file):
        if gene_id in annotations_data:
            merged_row = annotations_data[gene_id] + row
            writer.writerow(merged_row)

# 记录脚本结束时间
end_time = time.time()

# 计算脚本运行的总时间
elapsed_time = end_time - start_time
print(f"Merge completed, results saved to {output_file}")
print(f"Script execution time: {elapsed_time:.2f} seconds")
