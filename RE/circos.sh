
##  计算100K窗口内SNP数量
vcftools --gzvcf 186_filtered_vcftools.noContig.recode.vcf --SNPdensity 100000 --out 100K.SNP.density
#  After filtering, kept 12551267 out of a possible 12551267 Sites
#  Run Time = 41.00 seconds


# 计算100K窗口内INDEL数量
vcftools --gzvcf filtered_indel.vcf.gz --SNPdensity 100000 --out 100K.indel.density
#  After filtering, kept 22595531 out of a possible 22595531 Sites
#  Run Time = 719.00 seconds

# 生成绘图所需位点密度文件
grep -v "Contig" 100K.SNP.density.snpden  |  awk -F "\t" '{print $1,$2,($2+99999),$3}'  > 100K.SNP.noContig.density.txt
grep -v "Contig" 100K.indel.density.snpden  | awk -F "\t" '{print $1,$2,($2+99999),$3}'  > 100K.indel.noContig.density.txt

## 计算染色体长度
seqkit fx2tab -l -n -i Of.genome.fasta | grep "^Chr" |awk '{print $1"\t"$2}' > Of.genome.len

## 生成绘图所需染色体文件 7 列
awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}'  Of.genome.len > Of.karyotype.circos.txt

