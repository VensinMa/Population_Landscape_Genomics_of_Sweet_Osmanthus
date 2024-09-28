## 剔除多余的网脉木犀 最后保留了6个outgroup个体 O-WMMX O-DSMX O-MXL O-MZGH O-XYWJM O-ZNMX 和 194个内类群个体
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/
mkdir vcftools
vcftools --vcf all.filtered.snp.unanchor.vcf --keep keep_200sample.id --recode --recode-INFO-all --out ./vcftools/200samples_filtered.snp.unanchor

bcftools reheader -s 200sample_rename.id  ./vcftools/200samples_filtered.snp.unanchor.recode.vcf  -o ./vcftools/200samples_filtered_renamed.snp.unanchor.recode.vcf 


cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools

gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.recode.vcf" \
    --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered.vcf"
## [2024年9月24日 CST 01:32:37] org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration done. Elapsed time: 481.65 minutes.
## Runtime.totalMemory()=214748364800


gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered.vcf" \
    --exclude-filtered \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.gatk.vcf"

vcftools --vcf 200samples_filtered_renamed.snp.unanchor.final.gatk.vcf \
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
    --out 200samples_filtered_renamed.snp.unanchor.final.vcftools

vcftools --vcf 200samples_filtered_renamed.snp.unanchor.final.vcftools.recode.vcf \
    --max-missing 0.8 \
    --recode \
    --recode-INFO-all \
    --out 200samples_filtered_renamed.snp.unanchor.final.vcftools.nomissing


###  LD.prue
plink --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.recode.vcf  \
    --indep-pairwise 50 5 0.2 \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.recode.vcf \
    --allow-extra-chr \
    --set-missing-var-ids @:#
'''
11384901 variants and 200 people pass filters and QC.
Note: No phenotypes present.
Pruned 738495 variants from chromosome 27, leaving 109481.
Pruned 571440 variants from chromosome 28, leaving 85097.
Pruned 525337 variants from chromosome 29, leaving 78108.
Pruned 480921 variants from chromosome 30, leaving 71670.
Pruned 503252 variants from chromosome 31, leaving 75535.
Pruned 425157 variants from chromosome 32, leaving 66066.
Pruned 454359 variants from chromosome 33, leaving 70580.
Pruned 368305 variants from chromosome 34, leaving 57779.
Pruned 436567 variants from chromosome 35, leaving 67840.
Pruned 383854 variants from chromosome 36, leaving 60503.
Pruned 407635 variants from chromosome 37, leaving 64727.
Pruned 425639 variants from chromosome 38, leaving 62004.
Pruned 418648 variants from chromosome 39, leaving 64718.
Pruned 446051 variants from chromosome 40, leaving 67186.
Pruned 422543 variants from chromosome 41, leaving 65162.
Pruned 415039 variants from chromosome 42, leaving 64115.
Pruned 409556 variants from chromosome 43, leaving 65293.
Pruned 336186 variants from chromosome 44, leaving 54861.
Pruned 339825 variants from chromosome 45, leaving 51842.
Pruned 346383 variants from chromosome 46, leaving 53334.
Pruned 355035 variants from chromosome 47, leaving 56429.
Pruned 335785 variants from chromosome 48, leaving 54130.
Pruned 324289 variants from chromosome 49, leaving 48140.
Pruning complete.  9870301 of 11384901 variants removed.
Marker lists written to
/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.recode.vcf.prune.in
and
/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.recode.vcf.prune.out
.
'''

vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --recode-INFO-all --recode --max-missing 1 \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing

###  剔除外类群样本 NFM_1   O_DSMX  O_MXL   O_MZGH  O_XYWJM  只保留了浙南木犀
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.recode.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_5_outgroup.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/195samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing

vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/195samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.recode.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --minGQ 10 \
    --minQ 30 \
    --min-meanDP 6 \
    --max-missing 1 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/195samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.final


###  剔除外类群样本 NFM_1   O_DSMX  O_MXL   O_MZGH  O_XYWJM  只保留了浙南木犀
vcftools --vcf /root/workspace/njtree/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.recode.vcf \
    --recode --recode-INFO-all  --remove /root/workspace/njtree/removed_6_outgroup.id \
    --out /root/workspace/njtree/194samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing

vcftools --vcf /root/workspace/njtree/194samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.recode.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --minGQ 10 \
    --minQ 30 \
    --min-meanDP 6 \
    --max-missing 1 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out /root/workspace/njtree/194samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.final

    
