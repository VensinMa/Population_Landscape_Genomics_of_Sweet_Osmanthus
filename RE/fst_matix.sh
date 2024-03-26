#!/bin/bash

# 设置VCF文件和群体信息文件的路径
VCF_FILE="your_vcf_file.vcf"
POP_FILE="186sample.pop"

# 读取独特的群体名到数组
mapfile -t GROUPS < <(cut -f2 $POP_FILE | sort | uniq)

# 生成群体对组合
PAIRS=()
for i in "${!GROUPS[@]}"; do
    for j in $(seq $((i + 1)) ${#GROUPS[@]}); do
        PAIRS+=("${GROUPS[i]},${GROUPS[j]}")
    done
done

# 函数，用于计算两个群体间的FST并保存到文件
compute_fst() {
    IFS=',' read -r GROUP1 GROUP2 <<< "$1"
    TEMP_FILE="fst_${GROUP1}_${GROUP2}"
    vcftools --vcf "$VCF_FILE" \
        --weir-fst-pop <(awk -v grp="$GROUP1" '$2==grp {print $1}' $POP_FILE) \
        --weir-fst-pop <(awk -v grp="$GROUP2" '$2==grp {print $1}' $POP_FILE) \
        --out "$TEMP_FILE"

    # 提取FST值并追加到结果文件
    MEAN_FST=$(grep "Weir and Cockerham mean Fst estimate" "${TEMP_FILE}.log" | cut -d " " -f7)
    WEIGHTED_FST=$(grep "Weir and Cockerham weighted Fst estimate" "${TEMP_FILE}.log" | cut -d " " -f7)
    echo "$GROUP1,$GROUP2,$MEAN_FST" >> mean_fst_results.csv
    echo "$GROUP1,$GROUP2,$WEIGHTED_FST" >> weighted_fst_results.csv
    
    # 删除临时文件
    # rm "${TEMP_FILE}.log"
}

export -f compute_fst
export VCF_FILE POP_FILE

# 初始化结果文件并写入表头
echo "Group1,Group2,MeanFST" > mean_fst_results.csv
echo "Group1,Group2,WeightedFST" > weighted_fst_results.csv

# 使用parallel执行
parallel -j 16 compute_fst ::: "${PAIRS[@]}"

echo "FST计算完成，结果保存在mean_fst_results.csv和weighted_fst_results.csv"
