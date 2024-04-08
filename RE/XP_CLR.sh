cd /mnt/e/mwx/workspace/pop/xpclr

mkdir Chr01_23
######################### 提取单条染色体 ###################################
## 全部单条循环
seq -w 1 23 | parallel -j 30 vcftools --vcf ../186_filtered_vcftools.noContig.recode.vcf \
                                       --recode \
                                       --recode-INFO-all \
                                       --chr Chr{} \
                                       --out ./Chr01_23/186_filtered_vcftools.noContig.Chr{}
## 单条                                      
echo 12 | parallel -j 20 vcftools --vcf ../186_filtered_vcftools.noContig.recode.vcf \
                                  --recode \
                                  --recode-INFO-all \
                                  --chr Chr{} \
                                  --out ./Chr01_23/186_filtered_vcftools.noContig.Chr{}



#################################### 
for k in $(seq -w 1 23)
do
    awk -v chr="Chr$k" '$1==chr {print " "$1":"$2 "\t1\t" $2/100000000 "\t" $2 "\t" $4 "\t" $5 }' ./Chr01_23/186_filtered_vcftools.noContig.Chr${k}.recode.vcf > ./Chr01_23/Chr${k}.snp &
done


for k in  $(seq -w 1 23)
do xpclr --out Chr${k} \
         --format vcf \
         --input ./Chr01_23/186_filtered_vcftools.noContig.Chr${k}.recode.vcf \
         --samplesA west_samples.txt \
         --samplesB east_samples.txt \
         --map ./Chr01_23/Chr${k}.snp \
         --chr Chr${k} \
         --ld 0.95 \
         --maxsnps 200 \
         --size 2000 \
         --step 2000;
done

seq -w 1 23 | parallel -j [N] xpclr --out Chr{} \
                                 --format vcf \
                                 --input ./Chr01_23/186_filtered_vcftools.noContig.Chr{}.recode.vcf \
                                 --samplesA west_samples.txt \
                                 --samplesB east_samples.txt \
                                 --map ./Chr01_23/Chr{}.snp \
                                 --chr Chr{} \
                                 --ld 0.95 \
                                 --maxsnps 200 \
                                 --size 2000 \
                                 --step 2000

mkdir xpclr_res
seq -w 1 23 | parallel -j 30 xpclr --out xpclr_res/xpclr_Chr{}.out \
                                 --format vcf \
                                 --input ./Chr01_23/186_filtered_vcftools.noContig.Chr{}.recode.vcf \
                                 --samplesA west_samples.txt \
                                 --samplesB east_samples.txt \
                                 --chr Chr{} \
                                 --ld 0.95 \
                                 --maxsnps 200 \
                                 --size 2000 \
                                 --step 2000

186_filtered_vcftools.noContig.Chr01.recode.vcf
xpclr --out xpclr_Chr01 \
--format vcf \
--input ./Chr01_23/186_filtered_vcftools.noContig.Chr01.recode.vcf \
--samplesA west_samples.txt \
--samplesB east_samples.txt \
--chr Chr01 \
--ld 0.95 \
--maxsnps 200 \
--size 2000 \
--step 2000


                                 
