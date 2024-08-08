#!/bin/bash

# 输入目录和输出文件
input_dir="/mnt/e/mwx/workspace/186sample/population_freq_results"
output_file="combined_population_freq.csv"

# 初始化输出文件，添加表头行
echo -n "SNP" > $output_file

# 添加基因型列标题 (FREQ1/FREQ2, ALLELE1/ALLELE2)
first_file=$(ls ${input_dir}/*.frq | head -n 1)
genotype_line=$(awk 'NR==2 {print $5}' $first_file | sed 's/:/\//')
echo -n ",${genotype_line}" >> $output_file

# 添加群体名称（文件前缀）作为列标题
for file in ${input_dir}/*.frq; do
    group_name=$(basename $file .frq)
    echo -n ",${group_name}" >> $output_file
done
echo "" >> $output_file

# 处理每个SNP并将信息添加到输出文件
awk 'NR > 1' $first_file | while read line; do
    snp=$(echo $line | awk '{print $1 ":" $2}')
    echo -n "$snp" >> $output_file
    echo -n ",$genotype_line" >> $output_file

    for file in ${input_dir}/*.frq; do
        freq1=$(awk -v snp="$snp" '$1 ":" $2 == snp {print $6}' $file)
        echo -n ",$freq1" >> $output_file
    done
    echo "" >> $output_file
done
