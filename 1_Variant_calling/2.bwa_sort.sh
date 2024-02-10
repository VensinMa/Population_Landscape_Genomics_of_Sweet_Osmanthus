#!/bin/bash

# 设置原始数据目录
raw_data_dir="/public1/guop/mawx/workspace/wild_snpcalling/2.cleaned_data"

# 设置输出目录
output_dir="/public1/guop/mawx/workspace/wild_snpcalling/3.bwa_sam_bam"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir/sam" "$output_dir/bam" "$output_dir/sorted_bam"

# 设置日志文件
log_file="$output_dir/bwa_processing.log"

# 记录开始时间
echo "BWA processing started at $(date)" >> "$log_file"

# 构建参考基因组索引（如果尚未构建）
if [ ! -f "$reference_genome.bwt" ]; then
    echo "Building BWA index for reference genome..." >> "$log_file"
    bwa index "$reference_genome" >> "$log_file" 2>&1
fi

if [ ! -f "$reference_genome.fai" ]; then
    echo "Creating fasta index for reference genome..." >> "$log_file"
    samtools faidx "$reference_genome" >> "$log_file" 2>&1
fi

dict_file="${reference_genome%.*}.dict"
if [ ! -f "$dict_file" ]; then
    echo "Creating sequence dictionary for reference genome..." >> "$log_file"
    gatk CreateSequenceDictionary -R "$reference_genome" -O "$dict_file" >> "$log_file" 2>&1
fi

# 定义处理函数
process_sample() {
    fq1=$1
    fq2="${fq1/_1.clean.fq.gz/_2.clean.fq.gz}"  # 根据fq1构造fq2的文件名
    base_name=$(basename "$fq1" "_1.clean.fq.gz")  # 提取基本文件名，用于输出文件命名

    # 执行BWA-MEM比对
    echo "Starting BWA-MEM for $base_name" >> "$log_file"
    bwa mem -t 4 -M -R "@RG\tID:$base_name\tSM:$base_name\tLB:$base_name" "$reference_genome" "$fq1" "$fq2" > "$output_dir/sam/${base_name}.sam"
    echo "BWA-MEM completed for $base_name" >> "$log_file"

    # 转换SAM为BAM
    echo "Converting SAM to BAM for $base_name" >> "$log_file"
    samtools view -bS "$output_dir/sam/${base_name}.sam" -o "$output_dir/bam/${base_name}.bam"
    echo "SAM to BAM conversion completed for $base_name" >> "$log_file"

    # 排序BAM文件
    echo "Sorting BAM for $base_name" >> "$log_file"
    samtools sort "$output_dir/bam/${base_name}.bam" -o "$output_dir/sorted_bam/${base_name}.sorted.bam"
    echo "BAM sorting completed for $base_name" >> "$log_file"

    # 为排序后的BAM文件创建索引
    echo "Indexing sorted BAM for $base_name" >> "$log_file"
    samtools index "$output_dir/sorted_bam/${base_name}.sorted.bam"
    echo "Indexing completed for $base_name" >> "$log_file"

    echo "All processing steps completed for $base_name" >> "$log_file"
}

export -f process_sample
export raw_data_dir
export output_dir
export reference_genome
export log_file

# 导出必要的环境变量，以便parallel可以使用它们
find "$raw_data_dir" -name '*_1.clean.fq.gz' | sort | parallel -j 40 process_sample {}

echo "BWA processing completed at $(date)" >> "$log_file"
