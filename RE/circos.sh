cd /root/workspace/186sample/stat


##  计算100K窗口内SNP数量
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --SNPdensity 100000 --out 186_filtered.LD.pruned.noContig_100K.SNP.density
# After filtering, kept 1789805 out of a possible 1789805 Sites
# Run Time = 5.00 seconds



# 计算100K窗口内INDEL数量
vcftools --vcf /root/workspace/186sample/indel/186_filtered_indel_vcftools.noContig.recode.vcf --SNPdensity 100000 --out 186_filtered.LD.pruned.noContig_100K.indel.density
# After filtering, kept 979974 out of a possible 979974 Sites
# Run Time = 5.00 seconds


# 生成绘图所需位点密度文件
grep -v "Contig" 186_filtered.LD.pruned.noContig_100K.SNP.density.snpden  |  awk -F "\t" '{print $1,$2,($2+99999),$3}'  > 186_filtered.LD.pruned.noContig_100K.SNP.noContig.density.txt
grep -v "Contig" 186_filtered.LD.pruned.noContig_100K.indel.density.snpden  | awk -F "\t" '{print $1,$2,($2+99999),$3}'  > 186_filtered.LD.pruned.noContig_100K.indel.noContig.density.txt

## 计算染色体长度
seqkit fx2tab -l -n -i Of.genome.fasta | grep "^Chr" |awk '{print $1"\t"$2}' > Of.genome.len

## 生成绘图所需染色体文件 7 列
awk '{print "chr\t-\t"$1"\t"$1"\t0\t"$2"\tchr"NR}'  Of.genome.len > Of.karyotype.circos.txt

