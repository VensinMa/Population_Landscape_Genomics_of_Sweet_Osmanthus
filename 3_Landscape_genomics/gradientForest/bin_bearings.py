import pandas as pd
import numpy as np
import argparse

def main():
    # 设置命令行参数
    # 使用帮助 python bin_bearings.py input.csv output.csv
    parser = argparse.ArgumentParser(description='统计方位角每6度的分布数量')
    parser.add_argument('input_file', help='输入CSV文件路径')
    parser.add_argument('output_file', help='输出CSV文件路径')
    args = parser.parse_args()

    # 读取CSV文件
    df = pd.read_csv(args.input_file)

    # 确保方位角列存在
    if 'bearing' not in df.columns:
        raise ValueError("CSV文件中必须包含'bearing'列")

    # 将方位角范围从 [-180, 180] 转换为 [0, 360]
    df['bearing'] = np.where(df['bearing'] < 0, df['bearing'] + 360, df['bearing'])

    # 将方位角分为每6°一个区间
    bins = np.arange(0, 361, 6)  # 创建从0到360，每6度一个区间
    labels = bins[:-1]  # 使用每个区间的起始值作为标签
    df['bearing_bin'] = pd.cut(df['bearing'], bins=bins, labels=labels, include_lowest=True)

    # 统计每个分组的单元数
    bearing_counts = df['bearing_bin'].value_counts().sort_index()

    # 将结果保存到新的DataFrame
    result_df = pd.DataFrame({
        'bearing_bin': bearing_counts.index.astype(float),  # 转换为浮点型以便于后续处理
        'count': bearing_counts.values
    })

    # 保存结果到CSV文件
    result_df.to_csv(args.output_file, index=False)
    print(f"数据处理完成，结果已保存到 {args.output_file}")

if __name__ == "__main__":
    main()
