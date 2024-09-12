#!/bin/bash

# 设置输出目录
output_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/per_chromosome_vcf"
merged_vcf_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf"
log_dir="$merged_vcf_dir/logs"
mkdir -p "$merged_vcf_dir" "$log_dir" "$merged_vcf_dir/tmp"

# 记录脚本开始时间
echo "Merging VCF files started at $(date)" | tee -a "$log_dir/gatk_MergeVcfs_processing.log"

# 生成 VCF 文件列表 (每个文件都以 .vcf 为后缀)
vcf_list_file="$merged_vcf_dir/raw_vcf.list"
ls ${output_dir}/Superscaffold*.vcf > "$vcf_list_file"

# 确保 VCF 列表文件存在且非空
if [[ ! -s "$vcf_list_file" ]]; then
    echo "Error: VCF list file is empty or not found!" | tee -a "$log_dir/gatk_MergeVcfs_processing.log"
    exit 1
fi

# 设置合并后的输出 VCF 文件
merged_vcf="$merged_vcf_dir/all.merge_raw.vcf"

# 运行 GATK MergeVcfs
gatk --java-options "-Xmx10g -Djava.io.tmpdir=${merged_vcf_dir}/tmp" MergeVcfs \
    -I "$vcf_list_file" \
    -O "$merged_vcf" 2>&1 | tee -a "$log_dir/gatk_MergeVcfs.log"

# 检查合并结果
if [[ $? -eq 0 ]]; then
    echo "Merging VCF files completed successfully at $(date)" | tee -a "$log_dir/gatk_MergeVcfs_processing.log"
else
    echo "Error occurred during VCF merging at $(date)" | tee -a "$log_dir/gatk_MergeVcfs_processing.log"
    exit 1
fi
