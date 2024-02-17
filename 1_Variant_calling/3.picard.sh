#!/bin/bash

# 设置原始数据目录
raw_data_dir="/public1/guop/mawx/workspace/wild_snpcalling/2.cleaned_data"

# 设置输出目录
output_dir="/public1/guop/mawx/workspace/wild_snpcalling/3.bwa_sam_bam"

# 设置GATK输出目录
gatk_dir="/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf"

# 定义参考基因组文件
reference_genome="/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta"

# 创建输出目录（如果不存在）
mkdir -p "$output_dir/sam" "$output_dir/bam" "$output_dir/sorted_bam" "$output_dir/markdup" "$gatk_dir" "$output_dir/tmp"

# 设置日志文件
log_file="$output_dir/gatk_picard_HaplotypeCaller_processing.log"

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理函数
process_sample() {
    fq1=$1
    fq2="${fq1/_1.clean.fq.gz/_2.clean.fq.gz}"  # 根据fq1构造fq2的文件名
    base_name=$(basename "$fq1" "_1.clean.fq.gz")  # 提取基本文件名，用于输出文件命名

    # 记录样本处理开始时间
    echo "Processing $base_name started at $(date)" >> "$log_file"

    # Picard MarkDuplicates去重
    echo "Marking duplicates for $base_name at $(date)" >> "$log_file"
    java -jar /public1/guop/mawx/software/picard/picard.jar MarkDuplicates \
        -I "$output_dir/sorted_bam/${base_name}.sorted.bam" \
        -O "$output_dir/markdup/${base_name}.markdup.bam" \
        -M "$output_dir/markdup/${base_name}.markdup.metrics.txt" \
        -MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1000 \
        --REMOVE_DUPLICATES false \
        --TMP_DIR $output_dir/tmp
    echo "Marking duplicates completed for $base_name at $(date)" >> "$log_file"

    # samtools index对去重后的BAM文件建立索引
    echo "Indexing marked duplicates BAM for $base_name at $(date)" >> "$log_file"
    samtools index "$output_dir/markdup/${base_name}.markdup.bam"
    echo "Indexing marked duplicates BAM completed for $base_name at $(date)" >> "$log_file"

    # 运行GATK HaplotypeCaller进行变异检测
    echo "Running GATK HaplotypeCaller for $base_name at $(date)" >> "$log_file"
    gatk --java-options '-Xmx32g -XX:ParallelGCThreads=16 -Djava.io.tmpdir=./tmp' HaplotypeCaller \
        -R "$reference_genome" \
        -I "$output_dir/markdup/${base_name}.markdup.bam" \
        -O "$gatk_dir/${base_name}_raw.gvcf" \
        --native-pair-hmm-threads 16 \
        -ERC GVCF \
        --dont-use-soft-clipped-bases
    echo "GATK HaplotypeCaller completed for $base_name at $(date)" >> "$log_file"

    # 记录样本处理完成时间
    echo "Processing $base_name completed at $(date)" >> "$log_file"
}

export -f process_sample
export raw_data_dir
export output_dir
export gatk_dir
export reference_genome
export log_file

# 导出必要的环境变量，以便parallel可以使用它们
find "$raw_data_dir" -name '*_1.clean.fq.gz' | sort | parallel -j 40 process_sample {}

# 记录脚本完成时间
echo "Script completed at $(date)" >> "$log_file"