# 无外类群与含外类群 SNP基础过滤
vcftools --gzvcf  229_filtered_snp.vcf.gz  --min-alleles 2 --max-alleles 2 --minDP 5  --minGQ 10 --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05   --recode  --recode-INFO-all  --out 229_filtered_vcftools.vcf.gz


vcftools --gzvcf  224_filtered_snp.vcf.gz  --min-alleles 2 --max-alleles 2 --minDP 5  --minGQ 10 --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05   --recode  --recode-INFO-all  --out 224_filtered_vcftools.vcf.gz 
## After filtering, kept 13215558 out of a possible 172544728 Sites
## Run Time = 70892.00 seconds


# LD 过滤 生成保留和剔除位点的ID
plink --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --indep-pairwise 50 5 0.2 --out LD  --allow-extra-chr  --set-missing-var-ids @:#  
##  wc -l  LD.prune.in
## 1822433 LD.prune.in

# 去除 Contig位点 未定位到染色体上的位点
grep -v "Contig" LD.prune.in > LD.prune.noContig.in
##  wc -l  LD.prune.noContig.in
## 1821208 LD.prune.noContig.in 去除Contig位点后的数量
