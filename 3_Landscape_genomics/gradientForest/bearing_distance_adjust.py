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
    """处理距离单位转换（米→千米）并保留原始列"""
    if 'predDist' not in df.columns:
        raise ValueError("CSV文件中缺少必需的predDist列")
    df['predDist_adjusted'] = df['predDist'] / 1000  # 新增调整列
    return df

def main():
    parser = argparse.ArgumentParser(
        description='地理数据处理工具（方位角调整+距离单位转换）',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''\
        使用示例:
          基本用法:
            python bearing_distance_adjust.py input.csv output.csv --distance --bearing
            
          仅转换距离:
            python bearing_distance_adjust.py input.csv output_km.csv --distance
            
          仅调整方位角:
            python bearing_distance_adjust.py input.csv output_bearing.csv --bearing
            
          自动批量处理(Windows PowerShell):方法1
          
            # 在 PowerShell 中直接运行以下命令
            cmd /c "for %f in (""path\forwardOffset*.csv"") do python bearing_distance_adjust.py ""%f"" ""%~nf_adjusted.csv"" --distance --bearing "
            
          自动批量处理(Windows PowerShell):方法2
          
            Get-ChildItem "F:\offset数据\distance\csv\forwardOffset*.csv" | ForEach-Object {
            # 生成完整输出路径（与原文件同目录）
            $outputPath = Join-Path $_.Directory.FullName ($_.BaseName + "_adjusted.csv")
            # 执行处理命令
            python bearing_distance_adjust.py $_.FullName $outputPath --distance --bearing
            Write-Host "已处理: $($_.Name) → $($outputPath)"
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
        help='启用距离单位转换（新增predDist_adjusted列）'
    )
    
    parser.add_argument(
        '--bearing',
        action='store_true',
        help='启用方位角调整（新增bearing_adjusted列）'
    )

    args = parser.parse_args()

    try:
        if not args.distance and not args.bearing:
            raise ValueError("必须至少指定一个处理选项（--distance 或 --bearing）")

        df = pd.read_csv(args.input_file)
        print(f"成功读取: {args.input_file} ({len(df)} 条记录)")

        processed_ops = []
        if args.distance:
            df = process_distance(df)
            processed_ops.append('predDist_adjusted (米→千米)')
            
        if args.bearing:
            df = process_bearing(df)
            processed_ops.append('bearing_adjusted (方位角调整)')

        df.to_csv(args.output_file, index=False)
        print(f"\n处理完成: {args.output_file}")
        print("执行操作: " + ", ".join(processed_ops))
        print("当前列字段:")
        print(df.columns.tolist())

    except FileNotFoundError:
        print(f"错误：文件未找到 - {args.input_file}")
        sys.exit(1)
    except pd.errors.EmptyDataError:
        print("错误：空文件或无效格式")
        sys.exit(1)
    except Exception as e:
        print(f"运行错误: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
