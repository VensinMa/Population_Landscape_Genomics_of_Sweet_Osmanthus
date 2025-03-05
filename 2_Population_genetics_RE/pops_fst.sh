#!/bin/bash

source $(which env_parallel.bash)

# 输入文件路径
VCF_FILE="194samples_9674_adaptive_snps.recode.vcf"
POP_INFO_FILE="194sample.pop"

# 创建群体文件夹和日志文件夹
mkdir -p pop_files
mkdir -p log_files

## 第一步：根据总的群体信息文件创建每个群体的.pop文件 
awk '{print $1 > "pop_files/"$2".pop"}' "$POP_INFO_FILE"

# 生成组合
pop_files=(pop_files/*.pop)
combinations_file="combinations.txt"
> "$combinations_file"

for (( i=0; i<${#pop_files[@]}; i++ )); do
    for (( j=i+1; j<${#pop_files[@]}; j++ )); do
        echo "${pop_files[i]}:${pop_files[j]}" >> "$combinations_file"
    done
done

## 第二步：计算两两群体间FST
# 2.1：定义计算 FST 的函数
calculate_fst() {
    pair=$1
    IFS=":" read -r group1 group2 <<< "$pair"
    base_name1=$(basename "${group1%.*}")
    base_name2=$(basename "${group2%.*}")
    output_base="fst_${base_name1}_vs_${base_name2}"
    output_file="log_files/${output_base}"
    echo "Calculating FST for $base_name1 and $base_name2..."

    # 2.2：计算 FST
    vcftools --vcf "$VCF_FILE" --weir-fst-pop "$group1" --weir-fst-pop "$group2" --fst-window-size 100000 --fst-window-step 10000 --out "$output_file" >& "${output_file}.log"

    # 2.3：从日志文件提取 FST
    MEAN_FST=$(grep "Weir and Cockerham mean Fst estimate" "${output_file}.log" | awk '{print $NF}')
    WEIGHTED_FST=$(grep "Weir and Cockerham weighted Fst estimate" "${output_file}.log" | awk '{print $NF}')

    # 输出结果到 CSV 文件
    echo "${base_name1},${base_name2},${MEAN_FST},${WEIGHTED_FST}" >> fst_results.csv
}

export -f calculate_fst
export VCF_FILE

# 创建或清空结果文件，并添加标题行
echo "Group1,Group2,Mean_FST,Weighted_FST" > fst_results.csv

# 使用 GNU Parallel 并行计算 FST，从文件中读取组合 （注意调整线程数）
env_parallel -j 20 calculate_fst :::: "$combinations_file"

## 第三步：使用 R 脚本生成两个 FST 矩阵
Rscript -e "
generate_fst_matrices <- function(input_file, mean_output_file, weighted_output_file) {
  data <- read.csv(input_file, header=TRUE)
  unique_groups <- unique(c(data\$Group1, data\$Group2))
  
  # 初始化矩阵
  mean_matrix <- matrix(0, nrow=length(unique_groups), ncol=length(unique_groups), dimnames=list(unique_groups, unique_groups))
  weighted_matrix <- matrix(0, nrow=length(unique_groups), ncol=length(unique_groups), dimnames=list(unique_groups, unique_groups))

  # 填充矩阵
  for (i in 1:nrow(data)) {
    row_name <- as.character(data[i, 'Group1'])
    col_name <- as.character(data[i, 'Group2'])
    mean_matrix[row_name, col_name] <- as.numeric(data[i, 'Mean_FST'])
    mean_matrix[col_name, row_name] <- as.numeric(data[i, 'Mean_FST'])
    weighted_matrix[row_name, col_name] <- as.numeric(data[i, 'Weighted_FST'])
    weighted_matrix[col_name, row_name] <- as.numeric(data[i, 'Weighted_FST'])
  }

  # 写出矩阵
  write.csv(mean_matrix, mean_output_file, row.names=TRUE)
  write.csv(weighted_matrix, weighted_output_file, row.names=TRUE)
}

# 调用函数
generate_fst_matrices('fst_results.csv', 'mean_fst_matrix.csv', 'weighted_fst_matrix.csv')
"

echo "FST analysis and matrix generation completed successfully."
