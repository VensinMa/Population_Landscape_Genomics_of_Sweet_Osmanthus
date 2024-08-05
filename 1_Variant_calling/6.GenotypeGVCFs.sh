#!/bin/bash

# 设置combined_gvcf文件路径（来自上一步的输出）
combined_gvcf_path="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/combined_raw.gvcf"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 设置gatk文件目录
gatk_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf"

# 设置输出的VCF文件路径
output_vcf_path="$gatk_dir/combined_genotyped.vcf"

# 设置日志文件
log_file="$gatk_dir/gatk_GenotypeGVCFs_processing.log"

# 创建输出目录（如果不存在）
mkdir -p "$(dirname "$output_vcf_path")" "$gatk_dir/tmp"

# 记录脚本开始时间
echo "GenotypeGVCFs Script started at $(date)" | tee -a "$log_file"

# 运行GATK GenotypeGVCFs
echo "Running GATK GenotypeGVCFs at $(date)" | tee -a "$log_file"
gatk --java-options "-Xms200G -Xmx200G -XX:ParallelGCThreads=20 -Djava.io.tmpdir=${gatk_dir}/tmp" GenotypeGVCFs \
    -R "$reference_genome" \
    -V "$combined_gvcf_path" \
    -O "$output_vcf_path"  2>&1 | tee -a "$log_file"
echo "GATK GenotypeGVCFs completed at $(date)" | tee -a "$log_file"

# 记录脚本完成时间
echo "GenotypeGVCFs Script completed at $(date)" | tee -a "$log_file"
