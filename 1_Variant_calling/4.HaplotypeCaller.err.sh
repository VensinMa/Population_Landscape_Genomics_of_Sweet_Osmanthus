#!/bin/bash

# 设置输入markdup.bam文件目录
markdup_bam_dir="/public1/guop/mawx/workspace/wild_snpcalling/3.bwa_sam_bam/markdup"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 设置gatk文件目录
gatk_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf"

# 设置raw_gvcf文件输出目录
raw_gvcf_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/raw_gvcf"

# 创建输出目录（如果不存在）
mkdir -p "$gatk_dir" "$gatk_dir/tmp" "$raw_gvcf_dir"

# 设置日志文件
log_file="$gatk_dir/gatk_HaplotypeCaller_processing.log"

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理BAM文件的函数
process_HaplotypeCaller() {
    markdup_bam_path=$1
    base_name=$(basename "$markdup_bam_path" ".markdup.bam")  # 从BAM文件名提取基本文件名
    outdir=$(dirname "$markdup_bam_path")  # 获取markdup_bam_path的目录路径

    # 运行GATK HaplotypeCaller进行变异检测
    echo "Running GATK HaplotypeCaller for $base_name at $(date)" >> "$log_file"
    gatk --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${gatk_dir}/tmp" HaplotypeCaller \
        -R "$reference_genome" \
        -I "$markdup_bam_dir/${base_name}.markdup.bam" \
        -O "$raw_gvcf_dir/${base_name}_raw.gvcf" \
        --native-pair-hmm-threads 4 \
        -ERC GVCF
    echo "GATK HaplotypeCaller completed for $base_name at $(date)" >> "$log_file"
}

export -f process_HaplotypeCaller
export markdup_bam_dir
export reference_genome
export gatk_dir
export raw_gvcf_dir
export log_file

# 使用find命令和parallel并行处理markdup目录中的以下BAM文件
find "$markdup_bam_dir" -name '*.markdup.bam' | sort | grep -E "DST_2|SFZ_9|DYNC_6|DST_5|DYNC_2|DYNC_5|DYNC_3|DYNC_16|SFZ_7|DST_4|DST_9|DYNC_1|DST_8|DST_1|DST_7|DYNC_7|DYNC_10" | parallel -j 15 process_HaplotypeCaller {}

# 记录脚本完成时间
echo "GATK HaplotypeCaller Script completed at $(date)" >> "$log_file"

