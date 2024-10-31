#!/usr/bin/env python3

# 使用方法: python get_ancestral_allele.py input.vcf outgroupname1 outgroupname2

import sys

# 检查输入参数数量是否正确
if len(sys.argv) < 4:
    print("Usage: python get_ancestral_allele.py input.vcf outgroupname1 outgroupname2")
    sys.exit(1)

# 读取输入参数
input_vcf = sys.argv[1]
outgroup_name1 = sys.argv[2]
outgroup_name2 = sys.argv[3]

# 打开输入和输出文件
with open(input_vcf, "r") as vcf, \
     open(input_vcf.replace(".vcf", "_estsfs.txt"), "w") as estsfs_out, \
     open(input_vcf.replace(".vcf", "_estsfs.positions.txt"), "w") as positions_out:

    for line in vcf:
        if line.startswith("#CHROM"):
            samples = line.strip().split()
            # 获取外群样本的索引位置
            try:
                outgroup_index1 = samples.index(outgroup_name1)
                outgroup_index2 = samples.index(outgroup_name2)
            except ValueError:
                print("Error: One or both outgroup names not found in VCF header.")
                sys.exit(1)
            continue  # 继续读取下一个行

        elif line.startswith("#"):
            continue  # 忽略其他头部信息

        # 处理变异位点
        fields = line.strip().split("\t")
        CHROM, POS, ref_allele, alt_allele = fields[0], fields[1], fields[3], fields[4]

        # 仅处理单碱基替换（SNV），忽略InDel
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            allele_dic = {"A": 0, "C": 0, "G": 0, "T": 0}

            # 统计内群样本的等位基因频率
            ingroup_GT = "".join([i.split(":")[0] for i in fields[9:]])
            ref_count = ingroup_GT.count("0")
            alt_count = ingroup_GT.count("1")
            allele_dic[ref_allele] += ref_count
            allele_dic[alt_allele] += alt_count
            ingroup = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])

            # 重置等位基因计数字典用于外群1
            allele_dic = {key: 0 for key in allele_dic}
            outgroup1_genotype = fields[outgroup_index1].split(":")[0].split("/")
            if outgroup1_genotype[0] == outgroup1_genotype[1]:  # 排除杂合基因型
                if outgroup1_genotype[0] == "0":
                    allele_dic[ref_allele] += 2
                elif outgroup1_genotype[0] == "1":
                    allele_dic[alt_allele] += 2
            outgroup1 = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])

            # 重置等位基因计数字典用于外群2
            allele_dic = {key: 0 for key in allele_dic}
            outgroup2_genotype = fields[outgroup_index2].split(":")[0].split("/")
            if outgroup2_genotype[0] == outgroup2_genotype[1]:  # 排除杂合基因型
                if outgroup2_genotype[0] == "0":
                    allele_dic[ref_allele] += 2
                elif outgroup2_genotype[0] == "1":
                    allele_dic[alt_allele] += 2
            outgroup2 = ",".join([str(allele_dic["A"]), str(allele_dic["C"]), str(allele_dic["G"]), str(allele_dic["T"])])

            # 检查是否两个外群都没有有效数据
            if outgroup1 == "0,0,0,0" and outgroup2 == "0,0,0,0":
                continue  # 跳过该位点

            # 写入文件
            estsfs_out.write(f"{ingroup}\t{outgroup1}\t{outgroup2}\n")
            positions_out.write(f"{CHROM}\t{POS}\n")

print("祖先和派生等位基因信息已写入输出文件。")
