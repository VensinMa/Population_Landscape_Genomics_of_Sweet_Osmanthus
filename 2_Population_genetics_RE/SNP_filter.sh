## 剔除多余的网脉木犀 最后保留了6个outgroup个体 O-WMMX O-DSMX O-MXL O-MZGH O-XYWJM O-ZNMX 和 194个内类群个体
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/
vcftools --vcf all.filtered.snp.unanchor.vcf --keep keep_200sample.id --recode --recode-INFO-all --out ./vcftools/200samples_filtered.snp.unanchor

bcftools reheader -s 200sample_rename.id  ./vcftools/200samples_filtered.snp.unanchor.recode.vcf  -o ./vcftools/200samples_filtered_renamed.snp.unanchor.recode.vcf 


cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools
gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.recode.vcf" \
    --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered.vcf"

gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered.vcf" \
    --exclude-filtered \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.gatk.vcf"

vcftools --vcf 200samples_filtered_renamed.snp.unanchor.final.gatk.vcf  --min-alleles 2  --max-alleles 2  --minDP 5  --minGQ 10  --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05 --recode --recode-INFO-all --out 200samples_filtered_renamed.snp.unanchor.final.vcftools.recode.vcf




