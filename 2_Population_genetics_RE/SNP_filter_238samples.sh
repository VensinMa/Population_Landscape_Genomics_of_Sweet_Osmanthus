
## 保留所有个体  包括194个内类群个体
mkdir -p /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/238sample && cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/238sample

## vcftools --vcf all.filtered.snp.unanchor.vcf --keep keep_200sample.id --recode --recode-INFO-all --out ./vcftools/200samples_filtered.snp.unanchor
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.filtered.snp.unanchor.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --minGQ 10 \
    --minQ 30 \
    --min-meanDP 6 \
    --max-missing 0.8 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out 238samples_filtered.snp.unanchor.final.vcftools
