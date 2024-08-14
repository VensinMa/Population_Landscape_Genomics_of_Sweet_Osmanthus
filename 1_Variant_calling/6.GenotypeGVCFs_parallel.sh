#!/bin/bash

# 设置输入的 GVCF 文件路径
combined_gvcf_path="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/combined_raw.gvcf"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 设置输出目录
output_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/per_chromosome_vcf"
log_dir="$output_dir/logs"
mkdir -p "$output_dir" "$log_dir" "$output_dir/tmp"

# 记录脚本开始时间
echo "GenotypeGVCFs per chromosome script started at $(date)" | tee -a "$log_dir/gatk_GenotypeGVCFs_processing.log"

# 定义染色体名称列表 (根据你的染色体命名方式)
# chromosomes=("Superscaffold1" "Superscaffold2" "Superscaffold3" ... "Superscaffold23")
chromosomes=($(seq -f "Superscaffold%g" 1 23))

# 遍历每个染色体并行运行
for chrom in "${chromosomes[@]}"; do
    # 设置输出 VCF 文件路径
    output_vcf="${output_dir}/${chrom}.vcf"
    
    # 运行 GATK GenotypeGVCFs
    echo "Running GATK GenotypeGVCFs for $chrom at $(date)" | tee -a "$log_dir/gatk_GenotypeGVCFs_processing.log"
    gatk --java-options "-Xms50G -Xmx50G -XX:ParallelGCThreads=10 -Djava.io.tmpdir=${output_dir}/tmp" GenotypeGVCFs \
        -R "$reference_genome" \
        -V "$combined_gvcf_path" \
        -O "$output_vcf" \
        -L "$chrom" 2>&1 | tee -a "$log_dir/gatk_GenotypeGVCFs_${chrom}.log" &
    
    # 限制并发数，这里假设同时运行5个
    if [[ $(jobs -r -p | wc -l) -ge 25 ]]; then
        wait -n  # 等待任意一个任务完成
    fi
done

# 等待所有任务完成
wait

# 记录脚本完成时间
echo "GenotypeGVCFs per chromosome script completed at $(date)" | tee -a "$log_dir/gatk_GenotypeGVCFs_processing.log"

