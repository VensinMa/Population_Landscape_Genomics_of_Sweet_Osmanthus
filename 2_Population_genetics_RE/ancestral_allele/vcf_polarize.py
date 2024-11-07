#!/usr/bin/env python3
import sys
import pandas as pd

# 输入验证
if len(sys.argv) != 3:
    print("用法: python vcf_polarize.py <输入.vcf> <estsfs输出文件>")
    sys.exit(1)

if not sys.argv[1].endswith('.vcf'):
    print("错误: 输入文件必须是 .vcf 格式")
    sys.exit(1)

# 文件操作
try:
    in_vcf = open(sys.argv[1], "r")
    out_vcf = open(sys.argv[1].replace("vcf", "polarized.vcf"), "w")
    estsfs_out = open(sys.argv[2], "r")
except FileNotFoundError as e:
    print(f"错误: {e}")
    sys.exit(1)

# 读取祖先等位基因概率
ancestral_pval = []
for line in estsfs_out:
    if not line.startswith("0"):  # 跳过以0开头的行
        ancestral_pval.append(float(line.strip().split()[2]))  # 转换为浮点数以便比较
estsfs_out.close()

# 处理 VCF 文件逻辑
for i, line in enumerate(in_vcf):
    if line[0] == "#":
        out_vcf.write(line)  # 写入注释行
    else:
        x = line.strip().split("\t")
        CHROM = x[0]
        POS = x[1]
        ref_allele = x[3]
        alt_allele = x[4]
        
        # 统计等位基因计数
        allele_counts = {"A": 0, "C": 0, "G": 0, "T": 0}
        for genotype in x[9:]:
            gtype = genotype.split(":")[0]  # 获取基因型
            allele_counts[ref_allele] += gtype.count("0")
            allele_counts[alt_allele] += gtype.count("1")

        # 确定主要和次要等位基因
        major_allele = max(allele_counts, key=allele_counts.get)
        minor_allele = min(allele_counts, key=allele_counts.get)

        # 选择祖先等位基因
        ancestral_allele = major_allele if ancestral_pval[i] > 0.9 else minor_allele

        # 极性化
        if ancestral_allele == ref_allele:
            y = [gtype[0:3] for gtype in x[9:]]
        else:
            x[3], x[4] = ancestral_allele, ref_allele  # 交换参考和替代等位基因
            trans = str.maketrans("0/1-0/0-1/1", "1/0 1/1 0/0")
            y = [gtype[0:3].translate(trans) for gtype in x[9:]]

        # 写入输出文件
        out_vcf.write("\t".join(x[0:9]) + "\t" + "\t".join(y) + "\n")

# 关闭文件
in_vcf.close()
out_vcf.close()
