#!/usr/bin/env python3

# Usage: python vcf_to_estsfs.py input.vcf outgroupname1 [outgroupname2] [outgroupname3]

import sys

# 检查输入参数是否满足至少指定了1个外类群样本的要求
if len(sys.argv) < 3:
    print("Usage: python vcf_to_estsfs.py input.vcf outgroupname1 [outgroupname2] [outgroupname3]")
    sys.exit(1)

# 读取输入参数
input_vcf = sys.argv[1]
outgroup_names = sys.argv[2:]  # 支持1到3个外群名称
output_file_estsfs = input_vcf.replace(".vcf", "_estsfs_input.txt")
output_file_positions = input_vcf.replace(".vcf", "_estsfs.positions.txt")
log_file = input_vcf.replace(".vcf", "_processing.log")

# 初始化计数器
total_sites = 0
kept_sites = 0
skipped_sites = []  # 用于记录被跳过的位点信息

# 打开输入和输出文件
with open(input_vcf, "r") as vcf, \
     open(output_file_estsfs, "w") as estsfs_out, \
     open(output_file_positions, "w") as positions_out, \
     open(log_file, "w") as log_out:

    for line in vcf:
        if line.startswith("#"):  # 跳过注释行
            if line.startswith("#CHROM"):  # 找到样本头部
                samples = line.strip().split()
                # 获取所有外类群样本的索引位置
                outgroup_indices = []
                for name in outgroup_names:
                    try:
                        outgroup_indices.append(samples.index(name))
                    except ValueError:
                        print(f"Error: Outgroup name '{name}' not found in VCF header.")
                        sys.exit(1)
            continue  # 处理下一行

        total_sites += 1  # 每读取一行非注释行就增加总计数

        # 处理变异位点
        fields = line.strip().split("\t")
        CHROM, POS, ref_allele, alt_allele = fields[0], fields[1], fields[3], fields[4]
        site_info = f"{CHROM}:{POS}"

        # 仅处理二等位位点
        if len(ref_allele) == 1 and len(alt_allele) == 1:
            allele_dic = {"A": 0, "C": 0, "G": 0, "T": 0}

            # 跳过任何外类群样本中存在缺失基因型的位点
            if any(fields[idx].split(":")[0] in ["./.", "."] for idx in outgroup_indices):
                skipped_sites.append(f"Warning: Skipping site at {site_info} due to incomplete genotype format in outgroup sample.")
                continue

            # 获取外类群基因型，考虑基因型同时存在“/”和“|”两种分隔格式
            genotypes = [fields[idx].split(":")[0].replace("|", "/") for idx in outgroup_indices]

            # 跳过任何外类群样本中存在杂合基因型的位点 如果你只想保留纯合位点
            if any(gt[0] != gt[2] for gt in genotypes):
                skipped_sites.append(f"Warning: Skipping site at {site_info} due to heterozygous genotype in outgroup sample.")
                continue

            # 确保所有外类群样本的基因型一致 如果你只想保留在指定外类群基因型相同位点
            if len(set(genotypes)) > 1:
                skipped_sites.append(f"Warning: Skipping site at {site_info} due to inconsistent genotypes in outgroup samples.")
                continue

            # 统计内群样本的等位基因频率
            ingroup_GT = "".join([fields[i].split(":")[0] for i in range(9, len(fields)) if i not in outgroup_indices])
            allele_dic = {allele: 0 for allele in ["A", "C", "G", "T"]}  # 重置字典
            allele_dic[ref_allele] = ingroup_GT.count("0")
            allele_dic[alt_allele] = ingroup_GT.count("1")
            ingroup = ",".join(str(allele_dic[allele]) for allele in ["A", "C", "G", "T"])

            # 计算每个外类群样本的等位基因频率
            outgroups = []
            for genotype in genotypes:
                allele_dic = {allele: 0 for allele in allele_dic}  # 重置字典
                if genotype == "0/0":
                    allele_dic[ref_allele] = 1
                elif genotype == "1/1":
                    allele_dic[alt_allele] = 1
                else:
                    # 处理杂合基因型（如 "0/1" 或 "1/0"）
                    allele_dic[ref_allele] = 0
                    allele_dic[alt_allele] = 0
                outgroup_str = ",".join(str(allele_dic[allele]) for allele in ["A", "C", "G", "T"])
                outgroups.append(outgroup_str)

            # 写入文件，确保输出格式为“内群\t外群1 外群2 ... 外群N”
            estsfs_out.write(f"{ingroup}\t" + " ".join(outgroups) + "\n")
            positions_out.write(f"{CHROM}\t{POS}\n")
            kept_sites += 1  # 成功保留位点计数
        else:
            skipped_sites.append(f"Warning: Skipping site at {site_info} due to non-biallelic variant.")

    # 写入日志文件
    log_out.write("Skipped sites:\n")
    log_out.write(f"Total sites processed: {total_sites}\n")
    log_out.write(f"Successfully kept sites: {kept_sites}\n")

    for warning in skipped_sites:
        log_out.write(f"{warning}\n")

# 输出最终统计信息
print(f"VCF file has been converted to est-sfs input file.\nTotal sites processed: {total_sites}\nSuccessfully kept sites: {kept_sites}")
