#!/bin/bash

# 设置输入markdup.bam文件目录
markdup_bam_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/4.picard_markdup_bam/markdup"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 设置gatk文件目录
gatk_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/5.GATK_gvcf_vcf"

# 设置raw_gvcf文件输出目录
raw_gvcf_dir="${gatk_dir}/1.HaplotypeCaller_raw_gvcf"

# 创建输出目录（如果不存在）
mkdir -p "$gatk_dir/tmp" "$raw_gvcf_dir"

# 设置日志文件，添加时间戳（格式：YYYYMMDD_HHMMSS）
log_file="$gatk_dir/1.GATK_HaplotypeCaller_processing_$(date +"%Y%m%d_%H%M%S").log"

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理BAM文件的函数
process_HaplotypeCaller() {
    markdup_bam_path=$1
    base_name=$(basename "$markdup_bam_path" ".markdup.bam")  # 从BAM文件名提取基本文件名

    # 运行GATK HaplotypeCaller进行变异检测
    echo "Running GATK HaplotypeCaller for $base_name at $(date)" >> "$log_file"
    gatk --java-options "-Xms10G -Xmx20G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${gatk_dir}/tmp" HaplotypeCaller \
        -R "$reference_genome" \
        -I "$markdup_bam_path" \
        -O "$raw_gvcf_dir/${base_name}_raw.gvcf" \
        --native-pair-hmm-threads 4 \
        -ERC GVCF || { echo "Error in GATK HaplotypeCaller for $base_name" >> "$log_file"; return 1; }
    
    echo "GATK HaplotypeCaller completed for $base_name at $(date)" >> "$log_file"
}

export -f process_HaplotypeCaller
export markdup_bam_dir
export reference_genome
export gatk_dir
export raw_gvcf_dir
export log_file

# 使用find命令和parallel并行处理markdup目录中的所有BAM文件
find "$markdup_bam_dir" -name '*.markdup.bam' | sort | parallel -j 20 process_HaplotypeCaller {}

# 记录脚本完成时间
echo "GATK HaplotypeCaller Script completed at $(date)" >> "$log_file"
