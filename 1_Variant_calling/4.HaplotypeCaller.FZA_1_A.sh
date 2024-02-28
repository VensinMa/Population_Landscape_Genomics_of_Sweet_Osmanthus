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
mkdir -p "$gatk_dir" "$gatk_dir/tmp" "$raw_gvcf_dir" "$gatk_dir/status"

# 设置日志文件
log_file="$gatk_dir/gatk_HaplotypeCaller_processing.log"

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理BAM文件的函数
process_HaplotypeCaller() {
    local markdup_bam_path=$1
    local base_name=$(basename "$markdup_bam_path" ".markdup.bam")
    local status_file="$gatk_dir/status/${base_name}.status"  # 定义任务状态文件

    echo "Running GATK HaplotypeCaller for $base_name at $(date)" >> "$log_file"
    if gatk --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=4 -Djava.io.tmpdir=${gatk_dir}/tmp" HaplotypeCaller \
        -R "$reference_genome" \
        -I "$markdup_bam_path" \
        -O "$raw_gvcf_dir/${base_name}_raw.gvcf" \
        --native-pair-hmm-threads 4 \
        -ERC GVCF; then
        echo "SUCCESS" > "$status_file"
    else
        echo "FAILED" > "$status_file"
    fi
    echo "GATK HaplotypeCaller completed for $base_name at $(date)" >> "$log_file"
}

export -f process_HaplotypeCaller
export markdup_bam_dir
export reference_genome
export gatk_dir
export raw_gvcf_dir
export log_file

# 检查当前正在运行的相关进程数
current_jobs=$(ps -u $(whoami) -o comm= | grep 'java' | wc -l)

# 根据当前运行的进程数动态调整新的任务并行数
max_jobs=40
let "new_jobs=max_jobs-current_jobs"
if [ "$new_jobs" -le 0 ]; then
    echo "Already running $current_jobs jobs, no new jobs will be started." >> "$log_file"
else
    echo "Starting up to $new_jobs new jobs." >> "$log_file"
    find "$markdup_bam_dir" -name '*.markdup.bam' | sort | grep -A9999999 "FZA_1.markdup.bam" | parallel -j $new_jobs process_HaplotypeCaller {}
fi

# 定义重试失败任务的函数
retry_failed_tasks() {
    for status_file in $(find "$gatk_dir/status" -type f -name "*.status"); do
        if grep -q "FAILED" "$status_file"; then
            local base_name=$(basename "$status_file" ".status")
            local bam_path="$markdup_bam_dir/${base_name}.markdup.bam"
            echo "Retrying $base_name at $(date)" >> "$log_file"
            process_HaplotypeCaller "$bam_path"
        fi
    done
}

# 每小时检查并重试失败的任务
while true; do
    retry_failed_tasks
    sleep 3600  # 等待一小时
done

# 记录脚本完成时间
echo "GATK HaplotypeCaller Script completed at $(date)" >> "$log_file"
