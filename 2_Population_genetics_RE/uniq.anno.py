# -*- coding: utf-8 -*-

# This Python script processes an annotation file (e.g., LYG.genomic.pep.Nr5.annotations) by filtering out uninformative annotations associated with each transcript, based on specific exclusion keywords. 
# The goal is to retain only the most informative annotation for each transcript. Here's a breakdown of the script:

import re

# 定义需要排除的关键词
exclude_keywords = [
    "hypothetical protein", "uncharacterized protein",
    "Hypothetical predicted protein", "unnamed protein product"
]

def filter_annotations(file_path, output_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()

    result = []
    i = 0
    while i < len(lines):
        # 提取转录本ID
        transcript = lines[i].split()[0]

        # 获取当前转录本的注释行
        annotations = []
        while i < len(lines) and lines[i].split()[0] == transcript:
            annotations.append(lines[i])
            i += 1

        # 选择有效的注释
        selected_annotation = None
        for annotation in annotations:
            # 检查注释是否包含排除关键词
            if not any(keyword.lower() in annotation.lower() for keyword in exclude_keywords):
                selected_annotation = annotation
                break

        # 如果都包含排除关键词，保留第一条注释
        if selected_annotation is None:
            selected_annotation = annotations[0]

        # 将结果添加到最终列表
        result.append(selected_annotation)

    # 将结果写入输出文件
    with open(output_path, 'w') as f:
        f.writelines(result)

# 调用函数，传入输入文件路径和输出文件路径
input_file = 'LYG.genomic.pep.Nr5.annotations'  # 输入文件路径
output_file = 'LYG.genomic.pep.Nr5.uniq.annotations'  # 输出文件路径

filter_annotations(input_file, output_file)

