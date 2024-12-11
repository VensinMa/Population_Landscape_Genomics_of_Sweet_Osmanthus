## 剔除多余的网脉木犀 最后保留 194个桂花内类群个体
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/ && mkdir 194sample \
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/all.filtered.snp.unanchor.vcf \
    --keep /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/keep_194sample.id \
    --recode --recode-INFO-all \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.snp.unanchor

bcftools reheader -s /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample_rename.id  /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered.snp.unanchor.recode.vcf  -o /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered_renamed.snp.unanchor.recode.vcf


cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample  & mkdir tmp
### GATK 加标签
gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filtered_renamed.snp.unanchor.recode.vcf" \
    --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "SNP_filter" \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filter.vcf"
## 21:38:44.422 INFO  ProgressMeter - Traversal complete. Processed 141490137 total variants in 2690.0 minutes.
## 21:38:53.973 INFO  VariantFiltration - Shutting down engine
## [2024年10月1日 CST 21:38:53] org.broadinstitute.hellbender.tools.walkers.filters.VariantFiltration done. Elapsed time: 2,690.22 minutes.
#  Runtime.totalMemory()=214748364800

### GATK 过滤位点
gatk --java-options "-Xms200G -Xmx500G -XX:ParallelGCThreads=32 -Djava.io.tmpdir=./tmp" SelectVariants \
    -R "/public1/guop/mawx/workspace/wild_snpcalling/0.genome/LYG.hic.fasta" \
    -V "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_filter.vcf" \
    --exclude-filtered \
    -O "/public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_snp.nounanchor.renamed.filtered.vcf"

vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/194samples_snp.nounanchor.renamed.filtered.vcf \
    --min-alleles 2  --max-alleles 2  --minDP 5  --minGQ 10  --minQ 30 --min-meanDP 6 \
    --max-missing 0.8  --maf 0.05 --recode --recode-INFO-all \
    --out 194samples_snp.nounanchor.renamed.filtered.vcftools
## After filtering, kept 194 out of 194 Individuals
## Outputting VCF file...
## After filtering, kept 11452745 out of a possible 141490137 Sites
## Run Time = 30470.00 seconds

####################################  LD 过滤 --indep-pairwise 50 5 0.2  ########################################
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample
plink --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --indep-pairwise 50 5 0.2 \
    --out LD  \
    --allow-extra-chr  \
    --set-missing-var-ids @:#  
##  wc -l  LD.prune.in
##  LD.prune.in

sed 's/:/ /g' LD.prune.in > LD.prune.in.VCFTOOLS.txt

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions LD.prune.in.VCFTOOLS.txt \
    --recode --recode-INFO-all \
    --out 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned
##  After filtering, kept 194 out of 194 Individuals
##  Outputting VCF file...
##  After filtering, kept 1497179 out of a possible 11452745 Sites
##  Run Time = 3462.00 seconds


## 转bed与ped格式
plink --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
    --make-bed --allow-extra-chr \
    --out 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink

plink --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
      --recode --allow-extra-chr \
      --out 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink
plink --file 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink \
    --recodeA --allow-extra-chr  \
    --out 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink

## 使用plink将vcf转为raw
plink --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
      --recode --allow-extra-chr \
      --out  194samples_snp.nounanchor.renamed.filtered.vcftools.plink
plink --file 194samples_snp.nounanchor.renamed.filtered.vcftools.plink \
    --recodeA --allow-extra-chr  \
    --out 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.recodeA
sed '1d; s/NA/9/g' 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.recodeA.raw | \
awk '{$1=$2=$3=$4=$5=$6=""; print $0}' | \
sed 's/ \{6\}//g' > 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.lfmm

sed '1d; s/NA/9/g' 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.lfmm


####################################  提取基因型环境关联分析中的适应性位点  ########################################
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample

sed 's/:/ /g' 1766_adaptive_snps.id > 1766_adaptive_snps_vcftools.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf  \
    --positions 1766_adaptive_snps_vcftools.id \
    --recode --recode-INFO-all \
    --out 194samples_1766_adaptive_snps
##  After filtering, kept 194 out of 194 Individuals
##  Outputting VCF file...
##  After filtering, kept 1766 out of a possible 1497179 Sites
##  Run Time = 61.00 seconds






