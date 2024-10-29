#!/bin/bash

# 设置Picard软件的路径
picard_jar_path="/public1/guop/mawx/software/picard/picard.jar"

# 设置输入的sorted_bam文件夹路径
input_sorted_bam_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/3.bwa_sam_bam/sorted_bam"

# 设置输出文件夹路径
output_picard_dir="/public1/guop/mawx/workspace/cultivar_snpcalling/4.picard_markdup_bam"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 在output_picard_dir中创建必要的子目录
mkdir -p "$output_picard_dir/markdup" "$output_picard_dir/tmp" "$output_picard_dir/picard_metrics"

# 设置日志文件
log_file="$output_picard_dir/picard_markduplicates_processing_$(date +"%Y%m%d_%H%M%S").log"

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理BAM文件的函数
process_picard() {
    sorted_bam_path=$1
    base_name=$(basename "$sorted_bam_path" ".sorted.bam")  # 从BAM文件名提取基本文件名

    # 检查输入文件是否存在
    if [ ! -f "$sorted_bam_path" ]; then
        echo "Input file $sorted_bam_path does not exist. Skipping..." >> "$log_file"
        return 1
    fi

    # Picard MarkDuplicates去重
    echo "Marking duplicates for $base_name at $(date)" >> "$log_file"
    java -Xms20G -Xmx30G -jar "$picard_jar_path" MarkDuplicates \
        -I "$sorted_bam_path" \
        -O "$output_picard_dir/markdup/${base_name}.markdup.bam" \
        -M "$output_picard_dir/picard_metrics/${base_name}.metrics.txt" \
        --MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
        --REMOVE_DUPLICATES false \
        --ASSUME_SORTED true \
        --VALIDATION_STRINGENCY LENIENT \
        --TMP_DIR "$output_picard_dir/tmp" >> "$log_file" 2>&1 || { echo "MarkDuplicates failed for $base_name" >> "$log_file"; return 1; }
    echo "Marking duplicates completed for $base_name at $(date)" >> "$log_file"

    # samtools index对去重后的BAM文件建立索引
    echo "Indexing marked duplicates BAM for $base_name at $(date)" >> "$log_file"
    samtools index -@ 16 "$output_picard_dir/markdup/${base_name}.markdup.bam" >> "$log_file" 2>&1 || { echo "Indexing failed for $base_name" >> "$log_file"; return 1; }
    echo "Indexing marked duplicates BAM completed for $base_name at $(date)" >> "$log_file"
}

export -f process_picard
export picard_jar_path
export output_picard_dir
export reference_genome
export log_file

# 使用find命令和parallel并行处理输入目录中选择的BAM文件
# find "$input_sorted_bam_dir" -name '*.sorted.bam' |  sort | grep -E "D-HLJD-L|D-LZDG-L|D-MHDG-L|Y-CYYG-D|Y-SSYG-L|Y-WYYG-L|Y-XLYG-L|Y-XY-L|Y-XYSG-L|Y-XYWYG-L|Y-YLBZ-L|Y-YLL-L|Y-YLYS-L|Y-YMYYG-L|Y-YuLYS-L|Y-YX-L|Y-YYBYG-D|Y-YYG-S|Y-YZBZ-L|Y-YZ-L|Y-ZHLG-L|Y-ZiYG-L|Y-ZYG-L|ZS-L" | parallel -j 24 process_picard {}
find "$input_sorted_bam_dir" -name '*.sorted.bam' |  sort | grep -E "J-HCJG-H|J-LaYJG-L|Y-CG-L|Y-YLYS-L" | parallel -j 24 process_picard {}

# 记录脚本完成时间
echo "Script completed at $(date)" >> "$log_file"


