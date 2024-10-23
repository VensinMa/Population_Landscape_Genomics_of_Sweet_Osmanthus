cd /home/vensin/workspace/Annovar

awk '$1 == "intergenic" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic.snp.id
##  6075980 intergenic.snp.id


vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic.snp

plink --vcf 194samples_filtered.intergenic.snp.recode.vcf \
    --indep-pairwise 50 5 0.2 \
    --out LD  \
    --allow-extra-chr  \
    --set-missing-var-ids @:# 

sed 's/:/ /g' LD.prune.in > intergenic.LD.prune.in.snp.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic.LD.prune.in.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic.LD.prune

mkdir -p /home/vensin/workspace/intergenic_faststructure/result
cd /home/vensin/workspace/intergenic_faststructure

## vcf转bed格式
plink --vcf /home/vensin/workspace/Annovar/194samples_filtered.intergenic.LD.prune.recode.vcf \
    --make-bed   --keep-allele-order  --allow-extra-chr \
    --out /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune 

cd /home/vensin/workspace/intergenic_faststructure/result/
seq 2 10 | parallel -j 10 "structure.py -K {} \
    --input=/home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune  \
    --output=/home/vensin/workspace/intergenic_faststructure/result/LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > intergenic_LD_faststructure_K_{}.log 2>&1" &

chooseK.py --input intergenic_LD_faststructure_K

## Model complexity that maximizes marginal likelihood = 
## Model components used to explain structure in data = 
