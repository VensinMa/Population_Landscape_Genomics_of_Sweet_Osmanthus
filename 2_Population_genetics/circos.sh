
# 计算100K窗口内SNP数量
vcftools --gzvcf 224_filtered_snp.vcf.gz --SNPdensity 100000 --out 100K.SNP.density

# 生成绘图所需文件
awk -F "\t" '{print $1,$2,($2+99999),$3}' 100K.SNP.density.snpden  > 100K.SNP.density.txt
