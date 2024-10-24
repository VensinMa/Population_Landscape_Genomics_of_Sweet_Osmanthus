#!/bin/bash

# 设置原始数据目录
raw_data_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/1.raw_fastq"

# 设置输出目录
output_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/2.cleaned_data"  # 确保这是您的目标输出目录

# 设置质量报告输出目录
report_dir="$output_dir/reports"

# 检查并创建输出目录
mkdir -p "$output_dir" || { echo "Failed to create output directory: $output_dir"; exit 1; }

# 检查并创建质量报告的输出目录
mkdir -p "$report_dir" || { echo "Failed to create report directory: $report_dir"; exit 1; }

# 日志文件
log_file="$output_dir/fastp_processing.log"

# 记录开始时间
echo "Starting fastp processing at $(date)" | tee -a "$log_file"

# 定义处理函数
process_fastp() {
    fq1=$1
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"  # 对应的配对文件

    # 提取文件名前缀用于输出文件名
    filename=$(basename "$fq1")
    prefix=$(echo "$filename" | awk -F '_' '{print $1}')
    number=$(echo "$filename" | awk -F '_' '{print $2}')

    # 构建fastp命令
    fastp_cmd="fastp -w 4 -i $fq1 -I $fq2 -o $output_dir/${prefix}_${number}_1.clean.fq.gz -O $output_dir/${prefix}_${number}_2.clean.fq.gz --html $report_dir/${prefix}_${number}_fastp_report.html --json $report_dir/${prefix}_${number}_fastp_report.json"

    # 运行fastp命令
    echo "Processing $filename at $(date)" | tee -a "$log_file"
    eval "$fastp_cmd" &>> "$log_file"
}

export -f process_fastp  # 导出函数以供 parallel 使用
export output_dir report_dir log_file  # 导出变量以供 parallel 使用

# 使用find查找所有的fastq文件并使用parallel并行处理
find "$raw_data_dir" -name '*_1.fq.gz' | sort | parallel -j 20 process_fastp

# 等待所有并行任务完成
echo "All fastp processes completed at $(date)" | tee -a "$log_file"
