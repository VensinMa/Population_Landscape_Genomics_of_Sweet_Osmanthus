#!/bin/bash

# 输入文件路径
VCF_FILE="/root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf"
POP_INFO_FILE="/root/workspace/186sample/186sample.pop"
OUTPUT_MATRIX="fst_matrix.csv"

# 检查是否安装必要的工具
tools=(vcftools R awk)
for tool in "${tools[@]}"; do
  if ! command -v "$tool" &> /dev/null; then
    echo "Error: $tool is not installed. Please install it before running this script."
    exit 1
  fi
done

# 创建必要的文件夹
mkdir -p pop_files log_files

# 第一步：根据群体信息文件创建每个群体的.pop文件
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
    log_file="log_files/$output_base.log"

    # 计算FST
    if ! vcftools --vcf "$VCF_FILE" --weir-fst-pop "$group1" --weir-fst-pop "$group2" \
        --fst-window-size 100000 --fst-window-step 10000 --out "log_files/$output_base" &> "$log_file"; then
        echo "Error: FST calculation failed for $base_name1 and $base_name2" >> error.log
        return
    fi

    # 提取结果
    MEAN_FST=$(grep "Weir and Cockerham mean Fst estimate" "$log_file" | awk '{print $NF}')
    WEIGHTED_FST=$(grep "Weir and Cockerham weighted Fst estimate" "$log_file" | awk '{print $NF}')

    # 写入结果
    (
        flock -x 200
        echo "${base_name1},${base_name2},${MEAN_FST},${WEIGHTED_FST}" >> fst_results.csv
    ) 200>"fst_results.csv.lock"
}

# 创建结果文件
echo "Group1,Group2,Mean_FST,Weighted_FST" > fst_results.csv

# 并行计算FST
export -f calculate_fst
export VCF_FILE
env_parallel -j 20 calculate_fst ::: "${combinations[@]}"

# 第二步：生成FST矩阵
R_SCRIPT=$(cat <<EOF
# 读取 CSV 文件并生成 FST 矩阵
generate_fst_matrix <- function(input_file, output_file) {
  data <- read.csv(input_file, header=TRUE)
  unique_groups <- unique(c(data\$Group1, data\$Group2))
  matrix <- matrix(0, nrow=length(unique_groups), ncol=length(unique_groups), dimnames=list(unique_groups, unique_groups))
  for (i in 1:nrow(data)) {
    row_name <- as.character(data[i, "Group1"])
    col_name <- as.character(data[i, "Group2"])
    fst_value <- as.numeric(data[i, "Mean_FST"])
    matrix[row_name, col_name] <- fst_value
    matrix[col_name, row_name] <- fst_value
  }
  write.csv(matrix, output_file, row.names=TRUE)
  cat("FST matrix saved to", output_file, "\n")
}
# 执行函数
generate_fst_matrix("fst_results.csv", "$OUTPUT_MATRIX")
EOF
)

# 调用 R 脚本
Rscript -e "$R_SCRIPT"

# 脚本完成
echo "All processes are completed. FST matrix is saved in $OUTPUT_MATRIX."
