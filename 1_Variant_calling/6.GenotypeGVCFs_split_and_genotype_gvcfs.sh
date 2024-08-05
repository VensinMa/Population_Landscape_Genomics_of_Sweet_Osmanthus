#!/bin/bash

# 设置输入GVCF文件路径
input_gvcf="input.gvcf" # 根据你的实际GVCF文件路径调整

# 定义参考基因组文件
reference_genome="genome.final.fa"

# 设置输出目录
output_dir="/path/to/output_dir"
mkdir -p "$output_dir"

# 设置临时目录
tmp_dir="$output_dir/tmp"
mkdir -p "$tmp_dir"

# 设置日志文件
log_file="$output_dir/genotype_gvcfs_processing.log"

# 记录脚本开始时间
echo "GenotypeGVCFs Script started at $(date)" | tee -a "$log_file"

# 对每个染色体进行 GenotypeGVCFs 并行处理
for scaffold in Superscaffold{1..23}; do
    echo "Running GATK GenotypeGVCFs for $scaffold at $(date)" | tee -a "$log_file"
    gatk --java-options "-Xms200G -Xmx200G -XX:ParallelGCThreads=20 -Djava.io.tmpdir=$tmp_dir" GenotypeGVCFs \
        -L $scaffold \
        -R "$reference_genome" \
        -O "$output_dir/${scaffold}.raw.vcf.gz" \
        -V "$input_gvcf" 2>&1 | tee -a "$log_file" &
done

# 等待所有并行任务完成
wait

# 合并所有 VCF 文件
echo "Merging VCF files at $(date)" | tee -a "$log_file"
ls $output_dir/Superscaffold*.raw.vcf.gz > "$output_dir/all_genotype.list"
gatk --java-options "-Xms200G -Xmx200G -XX:ParallelGCThreads=20 -Djava.io.tmpdir=$tmp_dir" MergeVcfs \
    -I "$output_dir/all_genotype.list" \
    -O "$output_dir/final_combined.vcf.gz" 2>&1 | tee -a "$log_file"

# 记录脚本完成时间
echo "GenotypeGVCFs Script completed at $(date)" | tee -a "$log_file"
