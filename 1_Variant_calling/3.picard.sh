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

# 记录脚本开始时间
echo "Script started at $(date)" >> "$log_file"

# 定义处理BAM文件的函数
process_bam() {
    sorted_bam_path=$1
    base_name=$(basename "$sorted_bam_path" ".sorted.bam")  # 从BAM文件名提取基本文件名

    # Picard MarkDuplicates去重
    echo "Marking duplicates for $base_name at $(date)" >> "$log_file"
    java -jar /public1/guop/mawx/software/picard/picard.jar MarkDuplicates \
        -I "$sorted_bam_path" \
        -O "$output_dir/markdup/${base_name}.markdup.bam" \
        -M "$output_dir/markdup/${base_name}.markdup.metrics.txt" \
        --REMOVE_DUPLICATES false \
        --TMP_DIR $output_dir/tmp
    echo "Marking duplicates completed for $base_name at $(date)" >> "$log_file"

    # samtools index对去重后的BAM文件建立索引
    echo "Indexing marked duplicates BAM for $base_name at $(date)" >> "$log_file"
    samtools index "$output_dir/markdup/${base_name}.markdup.bam"
    echo "Indexing marked duplicates BAM completed for $base_name at $(date)" >> "$log_file"

    # 运行GATK HaplotypeCaller进行变异检测
    echo "Running GATK HaplotypeC(aller for $base_name at $(date)" >> "$log_file"
    gatk --java-options '-Xmx32g -XX:ParallelGCThreads=16 -Djava.io.tmpdir=./tmp' HaplotypeCaller \
        -R "$reference_genome" \
        -I "$output_dir/markdup/${base_name}.markdup.bam" \
        -O "$gatk_dir/${base_name}_raw.gvcf" \
        -ERC GVCF \
        --dont-use-soft-clipped-bases
    echo "GATK HaplotypeCaller completed for $base_name at $(date)" >> "$log_file"
}

export -f process_bam
export output_dir
export gatk_dir
export reference_genome
export log_file

# 使用parallel命令并行处理sorted_bam目录中的所有BAM文件
find "$output_dir/sorted_bam" -name '*.sorted.bam' | parallel -j 40 process_bam {}

# 记录脚本完成时间
echo "Script completed at $(date)" >> "$log_file"