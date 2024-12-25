#!/bin/bash

# 输入文件和目录
vcf_file="/home/vensin/workspace/Population_maf/194samples_9674_adaptive_snps.recode.vcf"
pop_file="/home/vensin/workspace/Population_maf/194sample.pop"

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"
combined_maf_file="combined_maf.csv"

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# 确保输入文件存在
if [ ! -f "$vcf_file" ]; then
    echo "VCF文件不存在：$vcf_file"
    exit 1
fi

if [ ! -f "$pop_file" ]; then
    echo "群体信息文件不存在：$pop_file"
    exit 1
fi

# 从群体信息文件中提取独特的群体，并保存到一个临时文件
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个独特的群体生成包含个体ID的文件，确保每行格式为 "个体ID 个体ID"
while read -r group; do
    grep "\b$group\b" "$pop_file" | awk '{print $1, $1}' > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# 初始化一个空的DataFrame用于存放结果
import pandas as pd
import os

# 获取所有.frq文件的列表
frq_files = [f for f in os.listdir(output_dir) if f.endswith('.frq')]

# 初始化一个空的DataFrame用于存放结果
combined_maf = pd.DataFrame()

# 初始化一个列表来存放所有文件中的SNP，用于保持顺序
all_snps = []

# 遍历所有群体文件并生成MAF
for group in $(cat "${temp_dir}/unique_groups.txt"); do
    echo "处理群体：$group"
    
    # 生成BED文件
    plink --vcf "$vcf_file" --keep "${temp_dir}/${group}.ids" --make-bed --allow-extra-chr --keep-allele-order --out "${output_dir}/${group}_maf" --set-missing-var-ids @:# 
    
    # 确认BED文件已生成
    if [ $? -eq 0 ]; then
        # 计算频率并保存为.frq文件
        plink --bfile "${output_dir}/${group}_maf" --allow-extra-chr --freq --out "${output_dir}/${group}_population_maf"
    else
        echo "生成BED文件失败：$group"
    fi

    # 读取.frq文件并保存MAF数据
    population_name="${group}"
    df = pd.read_csv(os.path.join(output_dir, f"{group}_population_maf.frq"), sep=r'\s+')
    df = df[['SNP', 'MAF']]
    df.rename(columns={'MAF': population_name}, inplace=True)

    # 如果是第一个文件，记录下SNP的顺序
    if not all_snps:
        all_snps = df['SNP'].tolist()

    # 合并数据
    if combined_maf.empty:
        combined_maf = df
    else:
        combined_maf = combined_maf.merge(df, on='SNP', how='outer')

# 对数据按SNP顺序进行排序
combined_maf['SNP'] = pd.Categorical(combined_maf['SNP'], categories=all_snps, ordered=True)
combined_maf = combined_maf.sort_values('SNP')

# 保存结果到CSV文件
combined_maf.to_csv("$combined_maf_file", index=False)
echo "合并的MAF数据已保存到 $combined_maf_file"
