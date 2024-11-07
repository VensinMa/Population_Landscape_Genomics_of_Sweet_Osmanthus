################################  NJ tree  #######################################

##
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools

python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf --fasta --output-prefix 200samples_filtered.snp.LD.pruned

setsid treebest nj -W -b 1000 200samples_filtered.snp.LD.pruned.min4.fasta > 200samples_filtered.snp.LD.pruned.treebest.out

sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 200samples_filtered.snp.LD.pruned.treebest.out | awk '{printf $0}' > 200samples_filtered.snp.LD.pruned.treebest.nwk

#################################  ML tree  #######################################

cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools
iqtree -s 200samples_filtered.snp.LD.pruned.min4.fasta -m MFP -B 1000 -m MFP -keep-ident  -T 20
#	-s					输入文件
#	-keep-ident			保留序列相同的样本
#	-st					指定输入序列的类别（当 iqtree 识别错误时使用）
#	-m					选择替换模型（substitution model）
#	-fast				加速建树
#	-nt					最大线程数

# 计算最佳模型并构树
iqtree -s 30252noouts.min4.fasta -m MFP -B 1000 -T 64
