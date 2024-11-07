#!/usr/bin/env python3

"""
indv_GT_stats.py
----------------
A script to calculate the genotype statistics (hom0, hom1, het, missing, total, and count_1) for each sample
from a VCF file. This script assumes the VCF contains genotype information for multiple individuals.

Usage:
    python indv_GT_stats.py <input_vcf_file> [<Type>]

Arguments:
    input_vcf_file    : The VCF file to analyze.
    Type              : (Optional) A string to filter variants by type, e.g., 'Radical', 'Moderately_Radical',
                        'Moderately_Conservative', 'Conservative'. If not provided, the script will process all variants.

Output:
    Prints the genotype statistics for each sample in the VCF file.

Example:
    python indv_GT_stats.py myfile.vcf
    python indv_GT_stats.py myfile.vcf Radical
"""

import sys
import vcf
from collections import defaultdict
import argparse

# Command-line argument parsing
def parse_args():
    parser = argparse.ArgumentParser(description="Calculate genotype statistics from a VCF file.")
    parser.add_argument("input_vcf", help="Input VCF file to analyze")
    parser.add_argument("Type", nargs="?", default="", help="Optional filter for variant type")
    return parser.parse_args()

# Main function
def main():
    # Parse command-line arguments
    args = parse_args()

    # Read VCF file using vcf.Reader
    vcf_file = args.input_vcf
    Type = args.Type  # Optional filter for variant type

    # Open the VCF file using the vcf module
    vcf_reader = vcf.Reader(open(vcf_file, "r"))
    sample_list = vcf_reader.samples

    # Initialize dictionaries for genotype counts
    GTdict = defaultdict(list)

    # Parse VCF file and store genotype data
    for record in vcf_reader:
        for sample in sample_list:
            GT = record.genotype(sample)["GT"]  # Get genotype for the sample
            GTdict[sample].append(GT)

    # Print header
    print("indv\thom0\thom1\thet\tmissing\ttotal\tcount_1")

    # Calculate statistics for each sample
    for sample in sample_list:
        hom0 = 0  # 0/0
        hom1 = 0  # 1/1
        het = 0   # 0/1 or 1/0
        missing = 0  # ./.
        count_1 = 0   # Count of 1 alleles

        # Iterate through genotypes to count occurrences
        for GT in GTdict[sample]:
            if GT == "0/0":
                hom0 += 1
            elif GT == "1/1":
                hom1 += 1
            elif GT in ("0/1", "1/0"):
                het += 1
                count_1 += 1  # Count 1 alleles in heterozygous
            elif GT == "./.":
                missing += 1
        
        total = hom0 + hom1 + het + missing
        count_1 += hom1 * 2  # Add 2 for each "1/1" (homozygous) genotype
        
        # Print the statistics for the current sample
        print(f"{sample}\t{hom0}\t{hom1}\t{het}\t{missing}\t{total}\t{count_1}")

if __name__ == "__main__":
    main()
