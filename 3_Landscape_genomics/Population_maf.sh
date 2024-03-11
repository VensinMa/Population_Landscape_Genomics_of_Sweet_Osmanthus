#!/bin/bash

# 输入文件: VCF文件和群体信息文件
vcf_file="224.filtered.LD.pruned.noContig.recode.vcf"
pop_info_file="224sample.pop" # 第一列为个体ID，第二列为群体名称

# 输出结果的文件夹
output_dir="population_maf_results"

# 如果输出文件夹不存在，则创建它
if [ ! -d "$output_dir" ]; then
    mkdir -p "$output_dir"
fi

# 为每个群体创建个体ID文件
awk '{print $1 > $2".ids"}' $pop_info_file

# 对每个群体计算次等位基因频率，并在后台执行
for pop in $(cut -f 2 -d ' ' $pop_info_file | sort | uniq)
do
    echo "Processing population: $pop"
    vcftools --vcf $vcf_file --keep ${pop}.ids --maf --out "${output_dir}/${pop}_maf" &
done

# 等待所有后台任务完成
wait

echo "All populations processed."




