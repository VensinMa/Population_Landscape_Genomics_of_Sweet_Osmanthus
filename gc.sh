# VCF文件转为BED BED转为RAW RAW转LFMM    0.01缺失过滤后 (6071370 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered_miss1_noscaffold.recode.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered_miss1_noscaffold.plink
  
plink --file /home/vensin/gc/311_filtered_miss1_noscaffold.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered_miss1_noscaffold.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered_miss1_noscaffold.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered_miss1_noscaffold.plink.recodeA.lfmm

awk '{print "Chr"$1":"$4}' 311_filtered_miss1_noscaffold.plink.map > 311_filtered_miss1_noscaffold.plink.ID


























# VCF文件转为BED BED转为RAW RAW转LFMM    0.01缺失过滤后 (6071370 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered_maf0.1_miss0.99.recode_noscaffold.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered_maf0.1_miss0.99.plink
  
plink --file /home/vensin/gc/311_filtered_maf0.1_miss0.99.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered_maf0.1_miss0.99.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered_maf0.1_miss0.99.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered_maf0.1_miss0.99.plink.recodeA.lfmm




# VCF文件转为BED BED转为RAW RAW转LFMM    LD过滤后 (2768836 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered.recode.LD.pruned.recode_noscaffold.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered.recode.LD.pruned.plink
  
plink --file /home/vensin/gc/311_filtered.recode.LD.pruned.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered.recode.LD.pruned.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered.recode.LD.pruned.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered.recode.LD.pruned.plink.recodeA.lfmm

# VCF文件转为BED BED转为RAW RAW转LFMM    零缺失过滤后 (974035 variants, 311 people)
plink --vcf /home/vensin/gc/311_filtered_maf0.1_nomissing_noscaffold.vcf\
      --recode --allow-extra-chr \
      --out  /home/vensin/gc/311_filtered_maf0.1_nomissing.plink
  
plink --file /home/vensin/gc/311_filtered_maf0.1_nomissing.plink \
    --recodeA --allow-extra-chr  \
    --out /home/vensin/gc/311_filtered_maf0.1_nomissing.plink.recodeA
    
sed '1d; s/NA/9/g' /home/vensin/gc/311_filtered_maf0.1_nomissing.plink.recodeA.raw | \
awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > /home/vensin/gc/311_filtered_maf0.1_nomissing.plink.recodeA.lfmm
