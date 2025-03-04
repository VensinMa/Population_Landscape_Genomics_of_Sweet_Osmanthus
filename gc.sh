# VCF文件转为BED BED转为RAW RAW转LFMM    0.01缺失过滤后 (8648736 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold.recode.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered.maf0.1.miss0.99.plink
  
plink --file /home/vensin/gc/311_filtered.maf0.1.miss0.99.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered.maf0.1.miss0.99.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered.maf0.1.miss0.99.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered.maf0.1.miss0.99.plink.recodeA.lfmm
awk '{print "Chr"$1":"$4}' 311_filtered.maf0.1.miss0.99.plink.map > 311_filtered.maf0.1.miss0.99.plink.ID


# VCF文件转为BED BED转为RAW RAW转LFMM    LD中性位点 (2768836 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.plink
  
plink --file /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.plink.recodeA.lfmm
awk '{print "Chr"$1":"$4}' 311_filtered.recode.LD.pruned.recode_noscaffold.plink.map > 311_filtered.recode.LD.pruned.recode_noscaffold.plink.ID



## 提取核心适应性位点
sed 's/:/ /g' /home/vensin/gc/lfmm_rda_intersection.csv > /home/vensin/gc/1002.coreSNP.id
vcftools --vcf /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold.recode.vcf \
    --positions 1002.coreSNP.id \
    --recode --recode-INFO-all \
    --out 311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP
##  After filtering, kept 1002 out of a possible 6071370 Sites
##  Run Time = 79.00 seconds


# VCF文件转为BED BED转为RAW RAW转LFMM    LD中性位点 (2768836 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.recode.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink
  
plink --file /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink.recodeA.lfmm
awk '{print "Chr"$1":"$4}' 311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink.map > 311_filtered.maf0.1.miss0.99_noscaffold_1002coreSNP.plink.ID




