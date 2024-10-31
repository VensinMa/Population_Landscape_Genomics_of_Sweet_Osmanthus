#!/usr/bin/env python3

# Usage: python get_ancestral_allele.py input.vcf outgroupname1 [outgroupname2] [outgroupname3]

import sys

# 检查输入参数是否满足至少指定了1个外类群样本的要求
if len(sys.argv) < 3:
    print("Usage: python get_ancestral_allele.py input.vcf outgroupname1 [outgroupname2] [outgroupname3]")
    sys.exit(1)

# 读取输入参数
input_vcf = sys.argv[1]
outgroup_names = sys.argv[2:]  # 支持1到3个外群名称
output_file_estsfs = input_vcf.replace(".vcf", "_estsfs_input.txt")
output_file_positions = input_vcf.replace(".vcf", "_estsfs.positions.txt")

# 打开输入和输出文件
with open(input_vcf, "r") as vcf, \
     open(output_file_estsfs, "w") as estsfs_out, \
     open(output_file_positions, "w") as positions_out:

    for line in vcf:
        if line.startswith("#CHROM"):
            samples = line.strip().split()
            # 获取所有外类群样本的索引位置
            outgroup_indices = []
            for name in outgroup_names:
                try:
                    outgroup_indices.append(samples.index(name))
                except ValueError:
                    print(f"Error: Outgroup name '{name}' not found in VCF header.")
                    sys.exit(1)
            continue  # 跳过表头行，处理下一个行

        elif line.startswith("#"):
            continue  # 忽略其他头部信息

        # 处理变异位点
        fields = line.strip().split("\t")
        CHROM, POS, ref_allele, alt_allele = fields[0], fields[1], fields[3], fields[4]

        # 仅处理二等位位点
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            allele_dic = {"A": 0, "C": 0, "G": 0, "T": 0}

            # 跳过任何外类群样本中存在缺失基因型的位点
            if any(fields[idx].split(":")[0] in ["./.", "."] for idx in outgroup_indices):
                continue

            # 获取外类群基因型，考虑基因型同时存在“/”和“|”两种分隔格式
            genotypes = [fields[idx].split(":")[0].replace("|", "/") for idx in outgroup_indices]

            # 跳过任何外类群样本中存在杂合基因型的位点 如果你只想保留纯合位点
            # if any(gt[0] != gt[2] for gt in genotypes):
            #     continue

            # 确保所有外类群样本的基因型一致 如果你只想保留在指定外类群基因型相同位点
            #if len(set(genotypes)) > 1:
            #    continue

            # 统计内群样本的等位基因频率
            ingroup_GT = "".join([i.split(":")[0] for i in fields[9:]])
            allele_dic = {allele: 0 for allele in allele_dic}  # 重置字典
            allele_dic[ref_allele] = ingroup_GT.count("0")
            allele_dic[alt_allele] = ingroup_GT.count("1")
            ingroup = ",".join(str(allele_dic[allele]) for allele in ["A", "C", "G", "T"])

            # 计算每个外类群样本的等位基因频率
            outgroups = []
            for genotype in genotypes:
                allele_dic = {allele: 0 for allele in allele_dic}  # 重置字典
                if genotype == "0/0":
                    allele_dic[ref_allele] = 2
                elif genotype == "1/1":
                    allele_dic[alt_allele] = 2
                else:
                    # 处理杂合基因型（如 "0/1" 或 "1/0"）
                    allele_dic[ref_allele] = 1
                    allele_dic[alt_allele] = 1
                outgroup_str = ",".join(str(allele_dic[allele]) for allele in ["A", "C", "G", "T"])
                outgroups.append(outgroup_str)

            # 如果存在缺失的外群基因型，则将缺失项填为 0,0,0,0
            # while len(outgroups) < 3:
            #     outgroups.append("0,0,0,0")

            # 写入文件，确保输出格式为“内群\t外群1 外群2 ... 外群N”
            estsfs_out.write(f"{ingroup}\t" + " ".join(outgroups) + "\n")
            positions_out.write(f"{CHROM}\t{POS}\n")

print("VCF file has been successfully converted to the est-sfs input file")
