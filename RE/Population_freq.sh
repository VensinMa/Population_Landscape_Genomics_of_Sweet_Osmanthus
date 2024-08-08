#!/bin/bash

# 输入文件和目录
vcf_file="186_filtered.LD.pruned.noContig.recode.vcf"
pop_file="186sample.pop"

# 输出目录
output_dir="population_freq_results"
temp_dir="${output_dir}/temp_ids"

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# 从群体信息文件中提取独特的群体，并保存到一个临时文件
cut -f2 $pop_file | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个独特的群体生成包含个体ID的文件
while read group; do
    grep "\b$group\b" $pop_file | awk '{print $1}' > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# 使用vcftools处理每个群体
while read group; do
    (
        echo "处理群体：$group"
        # 使用vcftools计算等位基因频率
        vcftools --vcf $vcf_file --keep "${temp_dir}/${group}.ids" --freq --out "${output_dir}/${group}_population_maf"
        if [ $? -eq 0 ]; then
            echo "完成群体：$group"
        else
            echo "计算等位基因频率失败：$group"
        fi
    ) &
done < "${temp_dir}/unique_groups.txt"

wait

echo "所有群体处理完毕。"
