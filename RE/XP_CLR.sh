cd /mnt/e/mwx/workspace/pop/xpclr

vcftools --vcf   ../186_filtered_vcftools.noContig.recode.vcf  --recode --recode-INFO-all --keep east.pop --out east.pop
vcftools --vcf   ../186_filtered_vcftools.noContig.recode.vcf  --recode --recode-INFO-all --keep west.pop --out west.pop

XPCLR -xpclr  -w1 0.0005 250 2500 $chr -p1 0.95

xpclr -xpclr  --format vcf west.pop.recode.vcf east.pop.recode.vcf   output.xpclr  -w1 0.005 200 2000 Chr01 -p0 0.95


for chr in {1..23}
do
  xpclr --format vcf west.pop.recode.vcf east.pop.recode.vcf output_chr${chr}.xpclr -w1 0.005 200 2000 ${chr} -p0 0.95
done

parallel -j 8 'xpclr -xpclr --format vcf --popA west.pop.recode.vcf --popB east.pop.recode.vcf --out output_Chr{}.xpclr --chr Chr{}  -w1 0.005 200 2000 -p0 0.95' ::: {01..23}
--chr 1 --maxsnps 600 --size 1000 --step 1000 --out

parallel -j 8 'python xpclr  --format vcf --popA west.pop.recode.vcf --popB east.pop.recode.vcf --out /root/workspace/186sample/xpclr/output_Chr{}.xpclr --chr Chr{}  --maxsnps 200 --size 2000 --step 2000 ' ::: {01..23}
--maxsnps 600 --size 1000 --step 1000 --out

parallel -j 25 'python /root/anaconda3/envs/xpclr/bin/xpclr  --format vcf --input ../186_filtered_vcftools.noContig.recode.vcf --samplesA  samplesA.id --samplesB  samplesB.id --out /root/workspace/186sample/xpclr/output_Chr{}.xpclr --chr Chr{}  --maxsnps 600 --size 1000 --step 1000 ' ::: {01..23}
