# -*- coding: utf-8 -*-
import pandas as pd
import argparse
import textwrap
import sys

def main():
    # 设置命令行参数
    parser = argparse.ArgumentParser(
        description='距离单位转换工具（米→千米）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        使用示例:
          基本用法:
            python convert_units.py input.csv output.csv
            
          自动批量处理(Windows):
            cmd /c "for %f in (""path\forwardOffset*.csv"") do python convert_units.py ""%f"" ""%~nf_km.csv"""
            
          处理含空格的文件:
            python convert_units.py "my input.csv" "my output.csv"
            
        输出文件格式:
          - 保留所有原始列
          - predDist列值会自动除以1000
          - 不添加额外索引列
        ''')
    )
    
    parser.add_argument(
        'input_file',
        metavar='INPUT_CSV',
        help='输入CSV文件路径（必须包含predDist列）'
    )
    
    parser.add_argument(
        'output_file',
        metavar='OUTPUT_CSV',
        help='输出CSV文件路径（建议使用_km后缀）'
    )

    args = parser.parse_args()

    try:
        # 读取CSV文件
        data = pd.read_csv(args.input_file)

        # 验证必要列存在
        if 'predDist' not in data.columns:
            raise ValueError("输入文件必须包含predDist列")

        # 执行单位转换
        data['predDist'] = data['predDist'] / 1000
        print(f"成功转换 {len(data)} 条记录的predDist值")

        # 保存结果
        data.to_csv(args.output_file, index=False)
        print(f"文件已保存至: {args.output_file}")

    except FileNotFoundError:
        print(f"错误：输入文件未找到 - {args.input_file}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("错误：输入文件为空或格式不正确")
        sys.exit(1)
    except Exception as e:
        print(f"处理过程中发生错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
