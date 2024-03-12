# 含外类群229 SNP基础过滤
vcftools --gzvcf  229_filtered_snp.vcf.gz  --min-alleles 2 --max-alleles 2 --minDP 5  --minGQ 10 --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05   --recode  --recode-INFO-all  --out 229_filtered_vcftools.vcf.gz
## After filtering, kept 12315778 out of a possible 183150846 Sites
## Run Time = 14132.00 seconds
sed '/^#CHROM/ s/_/-/g' 229_filtered_vcftools.vcf.gz.recode.vcf 

# 无外类群224 SNP基础过滤
vcftools --gzvcf  224_filtered_snp.vcf.gz  --min-alleles 2 --max-alleles 2 --minDP 5  --minGQ 10 --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05   --recode  --recode-INFO-all  --out 224_filtered_vcftools.vcf.gz 
## After filtering, kept 13215558 out of a possible 172544728 Sites
## Run Time = 70892.00 seconds

sed '/^#CHROM/ s/_/-/g' 229_filtered_vcftools.vcf.gz.recode.vcf > modified_229_filtered_vcftools.vcf.gz.recode.vcf

# LD 过滤 生成保留和剔除位点的ID
plink --vcf 229_filtered_vcftools.vcf.gz.recode.vcf --indep-pairwise 50 5 0.2 --out LD_229  --allow-extra-chr  --set-missing-var-ids @:#  
## wc -l  LD.prune.in
## 1728817

plink --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --indep-pairwise 50 5 0.2 --out LD  --allow-extra-chr  --set-missing-var-ids @:#  
##  wc -l  LD.prune.in
## 1822433
grep -v "Contig" LD.prune.in > LD.prune.noContig.in
## wc -l LD.prune.noContig.in
## 1727654 LD.prune.noContig.in
sed 's/:/ /g' LD.prune.noContig.in > vcftools_specific_noContig_position.txt

# 去除 Contig位点 未定位到染色体上的位点
grep -v "Contig" LD.prune.in > LD.prune.noContig.in
##  wc -l  LD.prune.noContig.in
## 1821208 LD.prune.noContig.in 去除Contig位点后的数量

## 根据LD.prune.noContig.in 要保留位点ID文件 提取 生成plink格式的vcf文件
plink --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --out LD.prune.noContig --extract LD.prune.in --recode vcf-iid  --allow-extra-chr  --keep-allele-order  --set-missing-var-ids @:#  
plink --vcf 229_filtered_vcftools.vcf.gz.recode.vcf --out LD.prune.noContig --extract LD.prune.in --recode vcf-iid  --allow-extra-chr  --keep-allele-order  --set-missing-var-ids @:#

# vcftools_specific_position.txt 使用vcftools 要保留位点ID文件提取vcftools格式的vcf文件 
vcftools --vcf 224_filtered_vcftools.vcf.gz.recode.vcf  --positions vcftools_specific_noContig_position.txt --recode --out 224_filtered.LD.pruned.noContig
##  保留下来的SNP数量 1821208 224_filtered.LD.pruned.noContig.recode.vcf

vcftools --vcf 229_filtered_vcftools.vcf.gz.recode.vcf  --positions vcftools_specific_noContig_position.txt --recode --out 229_filtered.LD.pruned.noContig










