cd /mnt/e/mwx/workspace/pop/genome

## 对参考基因组创建索引，生成五个文件（.amb .ann .bwt .pac .sa）
bwa index Osmanthus.genomic.fasta
### 生成fai索引
samtools faidx Osmanthus.genomic.fasta
### 生成dict索引
gatk CreateSequenceDictionary -R Osmanthus.genomic.fasta

 
cd /mnt/e/mwx/workspace/pop/
###  原始文件样品名称重命名
#根据输出文件格式要求，选择-Ov（不压缩输出）或者-Oz(bgzip压缩输出)参数
bcftools reheader -s rename.id 224_filtered_snp.vcf.gz -o 224_filtered_rename.vcf.gz 

###  剔除不确定样本
vcftools --gzvcf 224_filtered_rename.vcf.gz --recode --recode-INFO-all  --remove-indv remove.id --out 186_raw
#########################################  GATK 粗过滤  ######################################################
mkdir tmp
gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=16 -Djava.io.tmpdir=./tmp"  VariantFiltration  -R /mnt/e/mwx/workspace/pop/genome/Osmanthus.genomic.fasta -V 186_raw.recode.vcf --filter-expression "QUAL < 30.0 || QD < 2.0 || MQ < 40.0 || FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name 'gatk_filter' -O 186_filter.vcf.gz
gatk --java-options "-Xmx30g -Xms30g -XX:ParallelGCThreads=16 -Djava.io.tmpdir=./tmp"  SelectVariants  -R /mnt/e/mwx/workspace/pop/genome/Osmanthus.genomic.fasta -V 186_filter.vcf.gz --exclude-filtered  -O 186_filtered.vcf.gz

#########################################  vcftools 过滤  #####################################################

vcftools --gzvcf 224_filtered_rename.vcf.gz --remove remove.id --min-alleles 2  --max-alleles 2  --minDP 5  --minGQ 10  --minQ 30 --min-meanDP 6  --max-missing 0.8  --maf 0.05 --recode --recode-INFO-all --out 186_filtered_vcftools
## After filtering, kept 186 out of 224 Individuals
## Outputting VCF file...
## After filtering, kept 12557928 out of a possible 172544728 Sites
## Run Time = 15319.00 seconds


grep -v "^Contig" 186_filtered_vcftools.recode.vcf > 186_filtered_vcftools.noContig.recode.vcf
