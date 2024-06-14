#!/bin/bash

# 定义 fastp 报告文件所在目录
report_dir="/public1/guop/mawx/workspace/wild_snpcalling/2.cleaned_data/reports"

# 输出文件
output_file="$report_dir/fastp_summary.csv"

# 初始化输出文件并添加标题行
echo -e "Sample,Total Reads (Before),Total Bases (Before),Q20 Bases (Before),Q30 Bases (Before),GC Content (Before),Total Reads (After),Total Bases (After),Q20 Bases (After),Q30 Bases (After),GC Content (After),Reads Passed Filters,Reads with Low Quality,Reads with Too Many N,Reads Too Short" > "$output_file"

# 遍历所有 fastp 报告文件
for report in "$report_dir"/*.html; do
    # 提取样本名
    sample_name=$(basename "$report" "_fastp_report.html")
    
    # 提取 Before filtering 信息
    total_reads_before=$(grep "total reads:" "$report" | sed -n '1p' | awk -F '[<>]' '{print $5}')
    total_bases_before=$(grep "total bases:" "$report" | sed -n '1p' | awk -F '[<>]' '{print $5}')
    q20_bases_before=$(grep "Q20 bases:" "$report" | sed -n '1p' | awk -F '[<>]' '{print $5}' | awk '{print $1}')
    q30_bases_before=$(grep "Q30 bases:" "$report" | sed -n '1p' | awk -F '[<>]' '{print $5}' | awk '{print $1}')
    gc_content_before=$(grep "GC content:" "$report" | sed -n '1p' | awk -F '[<>]' '{print $5}')
    
    # 提取 After filtering 信息
    total_reads_after=$(grep "total reads:" "$report" | sed -n '2p' | awk -F '[<>]' '{print $5}')
    total_bases_after=$(grep "total bases:" "$report" | sed -n '2p' | awk -F '[<>]' '{print $5}')
    q20_bases_after=$(grep "Q20 bases:" "$report" | sed -n '2p' | awk -F '[<>]' '{print $5}' | awk '{print $1}')
    q30_bases_after=$(grep "Q30 bases:" "$report" | sed -n '2p' | awk -F '[<>]' '{print $5}' | awk '{print $1}')
    gc_content_after=$(grep "GC content:" "$report" | sed -n '2p' | awk -F '[<>]' '{print $5}')
    
    # 提取 Filtering result 信息
    reads_passed_filters=$(grep "reads passed filters:" "$report" | awk -F '[<>]' '{print $5}')
    reads_low_quality=$(grep "reads with low quality:" "$report" | awk -F '[<>]' '{print $5}')
    reads_too_many_N=$(grep "reads with too many N:" "$report" | awk -F '[<>]' '{print $5}')
    reads_too_short=$(grep "reads too short:" "$report" | awk -F '[<>]' '{print $5}')

    # 将结果写入输出文件
    echo -e "$sample_name,$total_reads_before,$total_bases_before,$q20_bases_before,$q30_bases_before,$gc_content_before,$total_reads_after,$total_bases_after,$q20_bases_after,$q30_bases_after,$gc_content_after,$reads_passed_filters,$reads_low_quality,$reads_too_many_N,$reads_too_short" >> "$output_file"
done

echo "Fastp report summary extraction completed. Results saved to $output_file."
