import os
import pandas as pd

# 输入目录和输出文件
input_dir = "/mnt/e/mwx/workspace/186sample/population_freq_results"
output_file = "combined_population_freq.csv"

# 获取所有 .frq 文件的列表
frq_files = [f for f in os.listdir(input_dir) if f.endswith('.frq')]

# 初始化字典来存储所有SNP信息
snp_data = {}

# 处理每个 .frq 文件
for frq_file in frq_files:
    file_path = os.path.join(input_dir, frq_file)
    
    # 从文件名获取群体名称，并去掉 "_population_maf" 部分
    group_name = os.path.splitext(frq_file)[0].replace('_population_maf', '')
    
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    # 遍历文件中的每一行数据
    for line in lines[1:]:  # 跳过标题行
        parts = line.strip().split()
        snp = f"{parts[0]}:{parts[1]}"
        allele_freqs = parts[4:]  # 提取所有等位基因频率对
        
        # 提取第一个等位基因频率
        first_allele, first_freq = allele_freqs[0].split(':')
        
        if snp not in snp_data:
            # 初始化FREQ列为 "allele1/allele2" 的形式
            snp_data[snp] = {"FREQ": f"{first_allele}/{allele_freqs[1].split(':')[0]}"}
        
        # 存储每个群体的第一个等位基因频率
        snp_data[snp][group_name] = first_freq

# 初始化 DataFrame 并写入表头
columns = ["SNP", "FREQ"] + [os.path.splitext(f)[0].replace('_population_maf', '') for f in frq_files]
df_list = []

# 填充 DataFrame
for snp, freqs in snp_data.items():
    row = {"SNP": snp, "FREQ": freqs["FREQ"]}
    for group_name in columns[2:]:
        row[group_name] = freqs.get(group_name, "NA")
    df_list.append(row)

df = pd.DataFrame(df_list, columns=columns)

# 保存结果到 CSV 文件
df.to_csv(output_file, index=False)

print(f"合并完成，结果保存在 {output_file}")
