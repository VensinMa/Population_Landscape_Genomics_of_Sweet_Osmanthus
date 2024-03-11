#!/bin/bash

# 输入的VCF文件
vcf_file="224.filtered.LD.pruned.noContig.recode.vcf"
# 群体信息文件
pop_info_file="224sample.pop"
# 输出目录
output_dir="output"

# 创建输出目录
mkdir -p "$output_dir"

# 读取群体信息并为每个群体创建个体ID列表文件
awk '{print $1 > "'$output_dir'/"$2".txt"}' $pop_info_file

# 定义一个函数来执行VCF提取和等位基因频率计算
process_population() {
    pop_name=$1
    pop_file="$output_dir/${pop_name}.txt"
    echo "Processing population: $pop_name"
    vcftools --vcf $vcf_file --keep $pop_file --freq --out "$output_dir/${pop_name}_af"
}

export -f process_population
export vcf_file
export output_dir

# 使用parallel并行执行process_population函数，限制同时运行的作业数为16
parallel -j 16 process_population ::: $(awk '{print $2}' $pop_info_file | sort | uniq)

echo "Parallel processing done. Now merging allele frequencies."

# 初始化合并文件
merged_file="$output_dir/merged_allele_frequencies.csv"
echo "Position" > $merged_file

# 为每个位点初始化行
declare -A pos_lines
for pop_file in $(ls $output_dir/*.frq); do
    pop_name=$(basename $pop_file _af.frq)
    while read -r line; do
        if [[ $line == CHROM* ]]; then
            continue
        fi
        pos=$(echo $line | awk '{print $1":"$2}')
        af=$(echo $line | awk '{print $5}')
        pos_lines[$pos]="${pos_lines[$pos]},$af"
    done < "$pop_file"
done

# 写入列名
pops=$(awk '{print $2}' $pop_info_file | sort | uniq | tr '\n' ',' | sed 's/,$//')
echo "Position,$pops" > $merged_file

# 将合并的行写入文件
for pos in "${!pos_lines[@]}"; do
    echo "$pos${pos_lines[$pos]}" >> $merged_file
done

echo "Merged allele frequencies are saved in $merged_file"
