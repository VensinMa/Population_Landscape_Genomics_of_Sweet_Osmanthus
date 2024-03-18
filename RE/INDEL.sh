cd /mnt/e/mwx/workspace/pop/genome

## 对参考基因组创建索引，生成五个文件（.amb .ann .bwt .pac .sa）
bwa index Osmanthus.genomic.fasta
### 生成fai索引
samtools faidx Osmanthus.genomic.fasta
### 生成dict索引
gatk CreateSequenceDictionary -R Osmanthus.genomic.fasta

 
cd /root/workspace/186sample/indel
###  原始文件样品名称重命名
#根据输出文件格式要求，选择-Ov（不压缩输出）或者-Oz(bgzip压缩输出)参数
bcftools reheader -s rename.id filtered_indel.vcf.gz -o 224_filtered_indel_rename.vcf.gz 

###  剔除不确定样本
vcftools --gzvcf 224_filtered_indel_rename.vcf.gz --recode --recode-INFO-all  --remove-indv remove.id --out 186_raw_indel

#########################################  vcftools 过滤  #####################################################

vcftools --gzvcf 224_filtered_indel_rename.vcf.gz  --remove remove.id --min-alleles 2  --max-alleles 2  --minDP 5  --minGQ 10  --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05 --recode --recode-INFO-all --out 186_filtered_indel_vcftools
##  After filtering, kept 186 out of 224 Individuals
##  Outputting VCF file...
##  After filtering, kept 980369 out of a possible 22595531 Sites
##  Run Time = 1101.00 seconds


grep -v "^Contig" 186_filtered_indel_vcftools.recode.vcf > 186_filtered_indel_vcftools.noContig.recode.vcf
grep -cv "#" 186_filtered_indel_vcftools.noContig.recode.vcf
## 979974








