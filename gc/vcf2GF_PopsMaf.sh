#!/bin/bash

# ===========================
# Shell脚本部分：处理群体MAF
# ===========================

# 输入文件和目录
vcf_file="/home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.recode.vcf"
pop_file="/home/vensin/gc/311sample.pop"

# 输出目录
output_dir="population_maf_results"
temp_dir="${output_dir}/temp_ids"

# 确保输出和临时目录存在
mkdir -p "$output_dir"
mkdir -p "$temp_dir"

# 从群体信息文件中提取独特的群体，并保存到一个临时文件
cut -f2 "$pop_file" | sort | uniq > "${temp_dir}/unique_groups.txt"

# 为每个独特的群体生成包含个体ID的文件，确保每行格式为 "个体ID 个体ID"
while read -r group; do
    grep "\b$group\b" "$pop_file" | awk '{print $1, $1}' > "${temp_dir}/${group}.ids"
done < "${temp_dir}/unique_groups.txt"

# 使用Plink处理每个群体
while read -r group; do
    (
        echo "=============================="
        echo "开始处理群体：$group"
        
        # 生成BED文件
        plink --vcf "$vcf_file" \
              --keep "${temp_dir}/${group}.ids" \
              --make-bed \
              --allow-extra-chr \
              --keep-allele-order \
              --out "${output_dir}/${group}_maf" \
              --set-missing-var-ids @:# 
        
        # 确认BED文件已生成
        if [ $? -eq 0 ]; then
            echo "BED文件生成成功：${output_dir}/${group}_maf.bed"
            
            # 计算频率
            plink --bfile "${output_dir}/${group}_maf" \
                  --allow-extra-chr \
                  --freq \
                  --out "${output_dir}/${group}_population_maf"
            
            if [ $? -eq 0 ]; then
                echo "频率计算成功：${output_dir}/${group}_population_maf.frq"
            else
                echo "频率计算失败：$group"
            fi
        else
            echo "生成BED文件失败：$group"
        fi
        
        echo "完成群体：$group"
        echo "=============================="
    ) &
done < "${temp_dir}/unique_groups.txt"

# 等待所有后台进程完成
wait

echo "所有群体处理完毕。"

# ====================================
# 嵌入第一个Python脚本：分析MAF结果
# ====================================

python3 << 'EOF1'
import os

def analyze_maf(output_dir):
    print("\n===== 开始分析MAF结果 =====")
    maf_files = [f for f in os.listdir(output_dir) if f.endswith("_population_maf.frq")]
    
    if not maf_files:
        print("未找到任何 .frq 文件。")
        return
    
    for maf_file in maf_files:
        file_path = os.path.join(output_dir, maf_file)
        if os.path.exists(file_path):
            with open(file_path, 'r') as f:
                lines = f.readlines()
                print(f"{maf_file} 包含 {len(lines)-1} 个SNP数据。")  # 减去标题行
        else:
            print(f"文件未找到：{maf_file}")

    print("===== MAF分析完成 =====\n")

def main():
    output_directory = "population_maf_results"
    analyze_maf(output_directory)

if __name__ == "__main__":
    main()
EOF1

# ====================================
# 嵌入第二个Python脚本：合并MAF结果
# ====================================

python3 << 'EOF2'
import pandas as pd
import os

def combine_maf(input_dir, output_file):
    print("===== 开始合并MAF结果 =====")
    frq_files = [f for f in os.listdir(input_dir) if f.endswith('.frq')]
    
    if not frq_files:
        print("未找到任何 .frq 文件用于合并。")
        return
    
    combined_maf = pd.DataFrame()
    all_snps = []
    
    for frq_file in frq_files:
        population_name = frq_file.split('_population_maf.frq')[0]
        file_path = os.path.join(input_dir, frq_file)
        
        try:
            df = pd.read_csv(file_path, sep=r'\s+')
            df = df[['SNP', 'MAF']]
            df.rename(columns={'MAF': population_name}, inplace=True)
            
            if not all_snps:
                all_snps = df['SNP'].tolist()
            
            if combined_maf.empty:
                combined_maf = df
            else:
                combined_maf = combined_maf.merge(df, on='SNP', how='outer')
            
            print(f"已处理文件：{frq_file}")
        except Exception as e:
            print(f"处理文件 {frq_file} 时出错：{e}")
    
    if not combined_maf.empty:
        combined_maf['SNP'] = pd.Categorical(combined_maf['SNP'], categories=all_snps, ordered=True)
        combined_maf.sort_values('SNP', inplace=True)
        combined_maf.set_index('SNP', inplace=True)
        combined_maf = combined_maf.T
        combined_maf.to_csv(output_file)
        print(f"===== 合并MAF结果已保存到 {output_file} =====\n")
    else:
        print("没有数据可供合并。")

def main():
    input_directory = "population_maf_results"
    output_csv = "GF_PopsMaf.csv"
    combine_maf(input_directory, output_csv)

if __name__ == "__main__":
    main()
EOF2

echo "所有操作已完成。"
