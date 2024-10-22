#!/bin/bash

# 设置参考基因组和输入/输出文件路径
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"
merged_vcf="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.merge_raw.vcf"
raw_snp_vcf="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.raw.snp.vcf"
filtered_snp_vcf="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.filter.snp.vcf"
final_snp_vcf="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.filtered.snp.vcf"
log_dir="./logs"
mkdir -p "$log_dir" "./tmp"

# 记录脚本开始时间
echo "SNP extraction, filtration, and final SNP extraction started at $(date)" | tee -a "$log_dir/gatk_snp_processing.log"

### 步骤 1：提取 SNP
echo "Step 1: Extracting SNP from $merged_vcf" | tee -a "$log_dir/gatk_snp_processing.log"
gatk --java-options "-Xms50G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "$reference_genome" \
    -V "$merged_vcf" \
    --select-type-to-include SNP \
    -O "$raw_snp_vcf" 2>&1 | tee -a "$log_dir/gatk_SelectVariants_snp.log"

# 检查提取结果
if [[ $? -eq 0 ]]; then
    echo "SNP extraction completed successfully." | tee -a "$log_dir/gatk_snp_processing.log"
else
    echo "Error during SNP extraction!" | tee -a "$log_dir/gatk_snp_processing.log"
    exit 1
fi

### 步骤 2：过滤 SNP（添加过滤标记）
echo "Step 2: Filtering SNP with quality filters" | tee -a "$log_dir/gatk_snp_processing.log"
gatk --java-options "-Xms50G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R "$reference_genome" \
    -V "$raw_snp_vcf" \
    --filter-expression "QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name 'SNP_filter' \
    -O "$filtered_snp_vcf" 2>&1 | tee -a "$log_dir/gatk_VariantFiltration_snp.log"

# 检查过滤结果
if [[ $? -eq 0 ]]; then
    echo "SNP filtration completed successfully." | tee -a "$log_dir/gatk_snp_processing.log"
else
    echo "Error during SNP filtration!" | tee -a "$log_dir/gatk_snp_processing.log"
    exit 1
fi

### 步骤 3：提取过滤通过的 SNP
echo "Step 3: Extracting filtered SNP" | tee -a "$log_dir/gatk_snp_processing.log"
gatk --java-options "-Xms50G -Xmx500G -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "$reference_genome" \
    -V "$filtered_snp_vcf" \
    --exclude-filtered \
    -O "$final_snp_vcf" 2>&1 | tee -a "$log_dir/gatk_SelectVariants_filtered_snp.log"

# 检查提取结果
if [[ $? -eq 0 ]]; then
    echo "Final SNP extraction completed successfully." | tee -a "$log_dir/gatk_snp_processing.log"
else
    echo "Error during filtered SNP extraction!" | tee -a "$log_dir/gatk_snp_processing.log"
    exit 1
fi

# 记录脚本结束时间
echo "SNP extraction, filtration, and final SNP extraction completed at $(date)" | tee -a "$log_dir/gatk_snp_processing.log"
