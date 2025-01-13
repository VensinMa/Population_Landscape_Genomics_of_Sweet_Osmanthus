#!/bin/bash

source $(which env_parallel.bash)

# 输入文件路径
VCF_FILE="194samples_9674_adaptive_snps.recode.vcf"
POP_INFO_FILE="194sample.pop"

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
    vcftools --vcf "$VCF_FILE" --weir-fst-pop "$group1" --weir-fst-pop "$group2" --fst-window-size 100000 --fst-window-step 10000 --out "log_files/$output_file" >& "log_files/$output_file.log"
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
env_parallel -j 8 calculate_fst ::: "${combinations[@]}"

# 使用R脚本生成FST矩阵
Rscript -e "
# 定义生成矩阵的函数
generate_fst_matrix <- function(input_file, output_file) {
  data <- read.csv(input_file, header=TRUE)
  unique_groups <- unique(c(data\$Group1, data\$Group2))
  matrix <- matrix(0, nrow=length(unique_groups), ncol=length(unique_groups), dimnames=list(unique_groups, unique_groups))
  for (i in 1:nrow(data)) {
    row_name <- as.character(data[i, 'Group1'])
    col_name <- as.character(data[i, 'Group2'])
    matrix[row_name, col_name] <- as.numeric(data[i, 'Mean_FST'])
    matrix[col_name, row_name] <- as.numeric(data[i, 'Mean_FST'])
  }
  write.csv(matrix, output_file)
}

# 调用函数
generate_fst_matrix('fst_results.csv', 'fst_matrix.csv')
"

echo "FST analysis and matrix generation completed successfully."
