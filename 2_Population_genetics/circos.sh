
##  计算100K窗口内SNP数量
vcftools --gzvcf 224_filtered_snp.vcf.gz --SNPdensity 100000 --out 100K.SNP.density
# After filtering, kept 172544728 out of a possible 172544728 Sites
# Run Time = 5679.00 seconds


# 计算100K窗口内INDEL数量

vcftools --gzvcf filtered_indel.vcf.gz --SNPdensity 100000 --out 100K.indel.density


# 生成绘图所需文件
awk -F "\t" '{print $1,$2,($2+99999),$3}' 100K.SNP.density.snpden  > 100K.SNP.density.txt


## 计算染色体长度

seqkit fx2tab -l -n -i Of.genome.fasta | grep "^Chr" |awk '{print $1"\t"$2}' > Of.genome.len

## 生成染色体文件 7 列

awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}'  Of.genome.len > Of.karyotype.circos.txt

