#!/bin/bash

# 设置输入raw_gvcf文件目录
raw_gvcf_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/raw_gvcf"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 设置gatk文件目录
gatk_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf"

# 设置combined_gvcf文件输出路径
combined_gvcf_path="$gatk_dir/combined_raw.gvcf"

# 设置日志文件
log_file="$gatk_dir/gatk_CombineGVCFs_processing.log"

# 创建输出目录（如果不存在）
mkdir -p "$combined_gvcf_path"

# 记录脚本开始时间
echo "CombineGVCFs Script started at $(date)" >> "$log_file"

# 获取所有GVCF文件路径
gvcf_files=$(find "$raw_gvcf_dir" -name '*.gvcf' | sort)

# 运行GATK CombineGVCFs
echo "Running GATK CombineGVCFs at $(date)" >> "$log_file"
gatk --java-options "-Xms200G -Xmx200G -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${gatk_dir}/tmp" CombineGVCFs \
    -R "$reference_genome" \
    --variant $gvcf_files \
    -O "$combined_gvcf_path"
echo "GATK CombineGVCFs completed at $(date)" >> "$log_file"

# 记录脚本完成时间
echo "CombineGVCFs Script completed at $(date)" >> "$log_file"

