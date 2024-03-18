#!/bin/bash

# 输入文件和目录
vcf_file="186_filtered.LD.pruned.noContig.recode.vcf"
pop_file="186sample.pop"

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# 从群体信息文件中提取独特的群体，并保存到一个临时文件
cut -f2 $pop_file | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个独特的群体生成包含个体ID的文件，确保每行格式为 "个体ID 个体ID"
while read group; do
    grep "\b$group\b" $pop_file | awk '{print $1, $1}' > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# 使用plink处理每个群体
while read group; do
    (
        echo "处理群体：$group"
        # 生成BED文件
        plink --vcf $vcf_file --keep "${temp_dir}/${group}.ids" --make-bed  --allow-extra-chr  --keep-allele-order  --out "${output_dir}/${group}_maf" --set-missing-var-ids @:# 
        # 确认BED文件已生成
        if [ $? -eq 0 ]; then
            # 计算频率
            plink --bfile "${output_dir}/${group}_maf" --freq --out "${output_dir}/${group}_population_maf"
        else
            echo "生成BED文件失败：$group"
        fi
        echo "完成群体：$group"
    ) &
done < "${temp_dir}/unique_groups.txt"

wait

echo "所有群体处理完毕。"
