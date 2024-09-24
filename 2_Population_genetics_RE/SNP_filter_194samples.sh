## 剔除多余的网脉木犀 最后保留 194个桂花内类群个体
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/ && mkdir 194sample
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.filtered.snp.unanchor.vcf --keep /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/keep_194sample.id --recode --recode-INFO-all --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.snp.unanchor

bcftools reheader -s /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample_rename.id  /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.snp.unanchor.recode.vcf  -o /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered_renamed.snp.unanchor.recode.vcf 


cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample

gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered_renamed.snp.unanchor.recode.vcf" \
    --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.vcf"
## [2024年9月24日 CST 01:32:37] org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration done. Elapsed time: 481.65 minutes.
## Runtime.totalMemory()=214748364800


gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.vcf" \
    --exclude-filtered \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered_renamed.snp.unanchor.final.gatk.vcf"

vcftools --vcf 200samples_filtered_renamed.snp.unanchor.final.gatk.vcf  --min-alleles 2  --max-alleles 2  --minDP 5  --minGQ 10  --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05 --recode --recode-INFO-all --out 200samples_filtered_renamed.snp.unanchor.final.vcftools.recode.vcf




