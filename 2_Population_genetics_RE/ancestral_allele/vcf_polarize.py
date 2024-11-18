#!/usr/bin/env python3
# usage: python vcf_polarize.py <input_vcf> <estsfs_output> <positions_file>

import sys
import argparse

# 创建解析器
def create_parser():
    parser = argparse.ArgumentParser(
        description="This script polarizes a VCF file based on ancestral allele probabilities from estsfs output."
    )
    parser.add_argument('input_vcf', type=str, help="Input VCF file to be polarized")
    parser.add_argument('estsfs_output', type=str, help="EstSFS output file containing ancestral allele probabilities")
    parser.add_argument('positions_file', type=str, help="File containing the positions of sites to be polarized from vcf_to_estsfs.py output")

    return parser

# 解析命令行参数
parser = create_parser()
args = parser.parse_args()

# 获取输入文件
try:
    in_vcf = open(args.input_vcf, "r")
    out_vcf = open(args.input_vcf.replace("vcf", "polarized.vcf"), "w")
    estsfs_out = open(args.estsfs_output, "r")
    positions_file = open(args.positions_file, "r")  # 读取 positions 文件
    ancestral_out = open(args.input_vcf.replace("vcf", "ancestral.txt"), "w")
except FileNotFoundError as e:
    print(f"Error: {e.strerror}: {e.filename}")
    sys.exit(1)

ancestral_out.write("CHROM\tPOS\tref_allele\tancestral_allele\n")

# 读取 positions 文件，获取需要处理的位点（染色体位置对）
positions = []
for line in positions_file:
    if line.strip():
        chr_pos = line.strip().split("\t")
        positions.append((chr_pos[0], int(chr_pos[1])))  # 存储为元组 (chromosome, position)
positions_file.close()

# 读取 estsfs 文件的祖先等位基因概率
ancestral_pval = []
for line in estsfs_out:
    if line[0] != "0":
        a = line.strip().split()
        ancestral_pval.append(a[2])  # 提取 P-major-ancestral 列
estsfs_out.close()

# 创建一个字典来记录 positions 中的位点位置及其索引
positions_dict = {f"{pos[0]}_{pos[1]}": idx for idx, pos in enumerate(positions)}

# 处理 VCF 文件，只处理 positions 文件中的位点
i = -1
for line in in_vcf:
    if line[0] == "#":
        out_vcf.write(line)  # 直接写入注释行
    else:
        i += 1
        x = line.strip("\n").split("\t")
        CHROM = x[0]
        POS = int(x[1])  # VCF 文件中的位点位置

        # 确保只处理 positions 文件中列出的位点
        if f"{CHROM}_{POS}" not in positions_dict:
            continue  # 如果当前位点不在 positions 文件中，跳过该位点

        # 获取该位点对应的 ancestral_pval 值
        estsfs_idx = positions_dict[f"{CHROM}_{POS}"]
        ancestral_probability = float(ancestral_pval[estsfs_idx])

        # 获取参考等位基因和变异等位基因
        ref_allele = x[3]
        alt_allele = x[4]
        allele_dic = {"A": 0, "G": 0, "C": 0, "T": 0}

        # 统计基因型信息
        for a in range(9, len(x)):
            ingroup_P = x[a]
            ingroup_genotype = ingroup_P[0:3].split("/")
            ingroup_0_count = ingroup_genotype.count("0")
            ingroup_1_count = ingroup_genotype.count("1")
            allele_dic[ref_allele] += ingroup_0_count
            allele_dic[alt_allele] += ingroup_1_count

        # 获取大等位基因和小等位基因
        major_allele = max(allele_dic, key=lambda k: allele_dic[k])
        allele_dic.pop(major_allele)
        minor_allele = max(allele_dic, key=lambda k: allele_dic[k])

        # 判断祖先等位基因
        if ancestral_probability > 0.9:
            ancestral_allele = major_allele
        else:
            ancestral_allele = minor_allele

        # 极性化：根据祖先等位基因调整参考和变异等位基因
        if ancestral_allele == ref_allele:
            y = [a[0:3] for a in x[9:]]  # 如果祖先等位基因与参考等位基因一致
        else:
            x[3] = ancestral_allele
            x[4] = ref_allele
            trans = str.maketrans("0/1-0/0-1/1", "1/0 1/1 0/0")  # 转换基因型
            y = [a[0:3].translate(trans) for a in x[9:]]

        # 写入极性化后的 VCF 文件
        out_vcf.write("\t".join(x[0:9]) + "\t" + "\t".join(y) + "\n")
        # 写入祖先等位基因信息
        ancestral_out.write(f"{CHROM}\t{POS}\t{ref_allele}\t{ancestral_allele}\n")

# 关闭文件
in_vcf.close()
out_vcf.close()
ancestral_out.close()
