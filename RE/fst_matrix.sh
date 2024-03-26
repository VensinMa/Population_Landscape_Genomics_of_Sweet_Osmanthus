#!/bin/bash

# 输入文件路径
VCF_FILE="/root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf"
POP_INFO_FILE="/root/workspace/186sample/186sample.pop"

# 创建群体文件夹和日志文件夹
mkdir -p pop_files
mkdir -p log_files

# 第一步：根据总的群体信息文件创建每个群体的.pop文件
awk '{print $1 > "pop_files/"$2".pop"}' $POP_INFO_FILE

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
    pair=($1)
    IFS=":" read -r -a groups <<< "$pair"
    group1="${groups[0]}"
    group2="${groups[1]}"
    base_name1=$(basename ${group1%.*})
    base_name2=$(basename ${group2%.*})
    output_base="fst_${base_name1}_vs_${base_name2}"
    output_file="log_files/${output_base}"
    echo "Calculating FST for $base_name1 and $base_name2..."
    
    # 第二步：计算FST
    vcftools --vcf $VCF_FILE --weir-fst-pop $group1 --weir-fst-pop $group2 --out $output_file
    OUTPUT=$(cat ${output_file}.log)

    # 第三步：提取FST估计值
    MEAN_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham mean Fst estimate" | awk '{print $7}')
    WEIGHTED_FST=$(echo "$OUTPUT" | grep "Weir and Cockerham weighted Fst estimate" | awk '{print $7}')
    
    # 输出结果到fst_results.txt，格式化为四列
    echo -e "${base_name1}\t${base_name2}\t$MEAN_FST\t$WEIGHTED_FST" >> fst_results.txt
}

export -f calculate_fst
export VCF_FILE

# 创建或清空结果文件，并添加标题行
echo -e "Group1\tGroup2\tMean_FST\tWeighted_FST" > fst_results.txt

# 使用GNU Parallel并行计算FST
parallel -j 4 calculate_fst ::: "${combinations[@]}"
