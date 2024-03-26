#!/bin/bash

source $(which env_parallel.bash)

# 输入文件路径
VCF_FILE="/root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf"
POP_INFO_FILE="/root/workspace/186sample/186sample.pop"

# 创建群体文件夹和日志文件夹
mkdir -p pop_files
mkdir -p log_files

# 第一步：根据总的群体信息文件创建每个群体的.pop文件
awk '{print $1 > "pop_files/"$2".pop"}' "$POP_INFO_FILE"

# 生成组合
pop_files=(pop_files/*.pop)
combinations=()

for (( i=0; i<${#pop_files[@]}; i++ )); do
    for (( j=i+1; j<${#pop_files[@]}; j++ )); do
        combinations+=("${pop_files[i]}:${pop_files[j]}")
    done
done

# 定义计算FST的函数
calculate_fst() {
    pair=$1
    IFS=":" read -r group1 group2 <<< "$pair"
    base_name1=$(basename "${group1%.*}")
    base_name2=$(basename "${group2%.*}")
    output_base="fst_${base_name1}_vs_${base_name2}"
    output_file="${output_base}"
    echo "Calculating FST for $base_name1 and $base_name2..."

    # 第二步：计算FST
    vcftools --vcf "$VCF_FILE" --weir-fst-pop "$group1" --weir-fst-pop "$group2" --fst-window-size 500000 --fst-window-step 50000 --out "log_files/$output_file" >& "log_files/$output_file.log"
    OUTPUT=$(cat "log_files/$output_file.log")

    # 第三步：从屏幕输出提取FST估计值
    MEAN_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham mean Fst estimate" | awk '{print $NF}')
    WEIGHTED_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham weighted Fst estimate" | awk '{print $NF}')

    # 输出结果到CSV文件，格式化为四列
    echo "${base_name1},${base_name2},${MEAN_FST},${WEIGHTED_FST}" >> fst_results.csv
}

export -f calculate_fst
export VCF_FILE

# 创建或清空结果文件，并添加标题行
echo "Group1,Group2,Mean_FST,Weighted_FST" > fst_results.csv

# 使用GNU Parallel并行计算FST
env_parallel -j 20 calculate_fst ::: "${combinations[@]}"
