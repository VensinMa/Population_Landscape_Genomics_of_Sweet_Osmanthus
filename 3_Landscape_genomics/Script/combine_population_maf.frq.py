import pandas as pd
import os

# 定义输入和输出路径
input_dir = 'population_maf_results'
output_file = 'combined_maf.csv'

# 获取所有.frq文件的列表
frq_files = [f for f in os.listdir(input_dir) if f.endswith('.frq')]

# 初始化一个空的DataFrame用于存放结果
combined_maf = pd.DataFrame()

# 遍历所有.frq文件
for frq_file in frq_files:
    # 解析群体名称，假设文件名格式为 "<群体名>_population_maf.frq"
    population_name = frq_file.split('_population_maf.frq')[0]
    
    # 读取.frq文件
    df = pd.read_csv(os.path.join(input_dir, frq_file), sep='\s+')
    
    # 仅保留 SNP 和 MAF 列
    df = df[['SNP', 'MAF']]
    
    # 重命名 MAF 列为群体名称
    df.rename(columns={'MAF': population_name}, inplace=True)
    
    # 合并数据
    if combined_maf.empty:
        combined_maf = df
    else:
        combined_maf = combined_maf.merge(df, on='SNP', how='outer')

# 将 SNP 列设置为索引
combined_maf.set_index('SNP', inplace=True)

# 转置 DataFrame，使得行为群体名称，列为 SNP
combined_maf = combined_maf.T

# 保存结果到 CSV 文件
combined_maf.to_csv(output_file)

