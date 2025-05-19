# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
import argparse
import textwrap
import sys

def process_bearing(df):
    """处理方位角转换 [-180, 180] → [0, 360]"""
    if 'bearing' not in df.columns:
        raise ValueError("CSV文件中缺少必需的bearing列")
    df['bearing_adjusted'] = np.where(df['bearing'] < 0, df['bearing'] + 360, df['bearing'])
    return df

def process_distance(df):
    """处理距离单位转换（米→千米）"""
    if 'predDist' not in df.columns:
        raise ValueError("CSV文件中缺少必需的predDist列")
    df['predDist'] = df['predDist'] / 1000
    return df

def main():
    parser = argparse.ArgumentParser(
        description='地理数据处理工具（方位角调整+距离单位转换）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        使用示例:
          基本用法:
            python geo_processor.py input.csv output.csv --distance --bearing
            
          仅转换距离:
            python geo_processor.py input.csv output_km.csv --distance
            
          仅调整方位角:
            python geo_processor.py input.csv output_bearing.csv --bearing
            
          自动批量处理(Windows PowerShell):
            Get-ChildItem *.csv | Foreach-Object {
              python geo_processor.py $_.FullName ($_.BaseName + "_processed.csv") --distance --bearing
            }
        ''')
    )
    
    parser.add_argument(
        'input_file',
        metavar='INPUT_CSV',
        help='输入CSV文件路径'
    )
    
    parser.add_argument(
        'output_file',
        metavar='OUTPUT_CSV',
        help='输出CSV文件路径'
    )
    
    parser.add_argument(
        '--distance',
        action='store_true',
        help='启用距离单位转换（米→千米）'
    )
    
    parser.add_argument(
        '--bearing',
        action='store_true',
        help='启用方位角调整（[-180,180]→[0,360]）'
    )

    args = parser.parse_args()

    try:
        # 参数验证
        if not args.distance and not args.bearing:
            raise ValueError("必须至少指定一个处理选项（--distance 或 --bearing）")

        # 读取数据
        df = pd.read_csv(args.input_file)
        original_rows = len(df)
        print(f"成功读取文件: {args.input_file} ({original_rows} 条记录)")

        # 执行处理流程
        processed_columns = []
        if args.distance:
            df = process_distance(df)
            processed_columns.append('predDist (米→千米)')
            
        if args.bearing:
            df = process_bearing(df)
            processed_columns.append('bearing (方位角调整)')

        # 保存结果
        df.to_csv(args.output_file, index=False)
        print(f"处理完成，已保存到: {args.output_file}")
        print(f"执行操作: {', '.join(processed_columns)}")
        print(f"输出文件列: {', '.join(df.columns)}")

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
