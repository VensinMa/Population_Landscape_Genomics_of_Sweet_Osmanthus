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






