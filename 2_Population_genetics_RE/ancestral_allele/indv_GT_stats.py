#!/usr/bin/env python3

"""
indv_GT_stats.py
----------------
A script to calculate the genotype statistics (hom0, hom1, het, missing, total, and count_1) for each sample
from a VCF file. This script assumes the VCF contains genotype information for multiple individuals.

Usage:
    python indv_GT_stats.py <input_vcf_file> [<Type>] [--output <output_file>]

Arguments:
    input_vcf_file    : The VCF file to analyze.
    Type              : (Optional) A string to filter variants by type, e.g., 'Radical', 'Moderately_Radical',
                        'Moderately_Conservative', 'Conservative'. If not provided, the script will process all variants.
    --output          : (Optional) Path to save the output file. If not provided, results are saved to 'indv_GT_stats_res.txt'.

Output:
    Saves the genotype statistics for each sample in the VCF file to the specified output file.

Example:
    python indv_GT_stats.py myfile.vcf
    python indv_GT_stats.py myfile.vcf Radical --output stats_output.txt
"""

import sys
import pandas as pd
import argparse
from collections import defaultdict

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate genotype statistics from a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file to analyze")
    parser.add_argument("Type", nargs="?", default="", help="Optional filter for variant type")
    parser.add_argument("--output", help="Optional output file to save results. Defaults to 'indv_GT_stats_res.txt'")
    return parser.parse_args()

# Main function
def main():
    # Parse command-line arguments
    args = parse_args()
    vcf_file = args.input_vcf
    variant_type = args.Type  # Optional filter for variant type
    output_file = args.output if args.output else "indv_GT_stats_res.txt"  # Default output file

    # Initialize dictionaries for genotype counts
    GTdict = defaultdict(list)
    sample_list = []

    # Read VCF file using pandas
    with open(vcf_file, "r") as file:
        lines = [line.strip() for line in file if not line.startswith("##")]
    header = lines[0].split("\t")
    df = pd.DataFrame([line.split("\t") for line in lines[1:]], columns=header)

    # Filter variants by type if specified
    if variant_type:
        df = df[df['INFO'].str.contains(variant_type, na=False)]

    # Extract sample names
    sample_list = header[9:]  # 从第10列开始就是样本数据

    # Parse genotype data for each sample
    for sample in sample_list:
        GTdict[sample] = df[sample].str.split(":", expand=True)[0].tolist()  # 仅提取基因型（GT字段）

    # Prepare output (either to specified file or default file)
    output_lines = []
    output_lines.append("indv\thom0\thom1\thet\tmissing\ttotal\tcount_1")

    # Calculate statistics for each sample
    for sample in sample_list:
        hom0 = 0  # 0/0
        hom1 = 0  # 1/1
        het = 0   # 0/1 or 1/0
        missing = 0  # ./.
        count_1 = 0   # Count of 1 alleles

        # Count genotype occurrences
        for GT in GTdict[sample]:
            if GT == "0/0":
                hom0 += 1
            elif GT == "1/1":
                hom1 += 1
            elif GT in ("0/1", "1/0"):
                het += 1
                count_1 += 1  # Count 1 alleles in heterozygous genotypes
            elif GT == "./.":
                missing += 1
        
        total = hom0 + hom1 + het + missing
        count_1 += hom1 * 2  # Each "1/1" genotype contributes 2 to count_1
        
        # Append statistics for the sample to output_lines
        output_lines.append(f"{sample}\t{hom0}\t{hom1}\t{het}\t{missing}\t{total}\t{count_1}")

    # Write results to the output file
    with open(output_file, "w") as f:
        f.write("\n".join(output_lines) + "\n")
    print(f"Results saved to '{output_file}'")

if __name__ == "__main__":
    main()
