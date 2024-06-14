#!/bin/bash

# 设置原始数据目录
raw_data_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/fastq.gz/raw_data"

# 设置输出目录
clean_data_dir="/public1/guop/mawx/workspace/wild_snpcalling/2.cleaned_data"  # 确保这是您的目标输出目录

# 设置质量报告输出目录
report_dir="$clean_data_dir/reports"

# 检查并创建输出目录
mkdir -p "$clean_data_dir" || { echo "Failed to create output directory: $clean_data_dir"; exit 1; }

# 检查并创建质量报告的输出目录
mkdir -p "$report_dir" || { echo "Failed to create report directory: $report_dir"; exit 1; }

# 日志文件
log_file="$clean_data_dir/fastp_processing.log"

# 记录开始时间
echo "Starting fastp processing at $(date)" | tee -a "$log_file"

# 定义处理函数
process_sample() {
    fq1=$1
    fq2="${fq1/_1.fq.gz/_2.fq.gz}"  # 对应的配对文件

    # 提取文件名前缀用于输出文件名
    filename=$(basename "$fq1")
    prefix=$(echo "$filename" | awk -F '_' '{print $1}')
    number=$(echo "$filename" | awk -F '_' '{print $2}')

    # 构建fastp命令
    fastp_cmd="fastp -w 4 -i $fq1 -I $fq2 -o $clean_data_dir/${prefix}_${number}_1.clean.fq.gz -O $clean_data_dir/${prefix}_${number}_2.clean.fq.gz --html $report_dir/${prefix}_${number}_fastp_report.html --json $report_dir/${prefix}_${number}_fastp_report.json"

    # 使用后台进程运行fastp命令，并行处理
    echo "Processing $filename at $(date)" | tee -a "$log_file"
    eval "$fastp_cmd" &>> "$log_file"
}

export -f process_sample
export output_dir
export report_dir
export log_file

# 导出必要的环境变量，以便parallel可以使用它们
find "$raw_data_dir" -name '*_1.fq.gz' | sort | parallel -j 40 process_sample {}

# 记录结束时间
echo "All fastp processes completed at $(date)" | tee -a "$log_file"
