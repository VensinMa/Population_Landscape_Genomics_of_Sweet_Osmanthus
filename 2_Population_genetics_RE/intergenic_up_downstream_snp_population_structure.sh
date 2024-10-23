cd /home/vensin/workspace/Annovar

##  awk '$1 == "intergenic" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic.snp.id
##  6075980 intergenic.snp.id

awk '$1 == "intergenic" || $1 == "upstream" || $1 == "downstream" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic_up_downstream.snp.id
##  8073005 intergenic_up_downstream.snp.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic_up_downstream.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic_up_downstream.snp

plink --vcf 194samples_filtered.intergenic_up_downstream.snp.recode.vcf \
    --indep-pairwise 50 5 0.2 \
    --out LD  \
    --allow-extra-chr  \
    --set-missing-var-ids @:# 

sed 's/:/ /g' LD.prune.in > intergenic_up_downstream.LD.prune.in.snp.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic_up_downstream.LD.prune.in.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic_up_downstream.LD.prune

mkdir -p /home/vensin/workspace/intergenic_up_downstream_faststructure/result
cd /home/vensin/workspace/intergenic_up_downstream_faststructure

## vcf转bed格式
plink --vcf /home/vensin/workspace/Annovar/194samples_filtered.intergenic_up_downstream.LD.prune.recode.vcf \
    --make-bed   --keep-allele-order  --allow-extra-chr \
    --out /home/vensin/workspace/intergenic_up_downstream_faststructure/194samples_filtered.intergenic_up_downstream.LD.prune 

cd /home/vensin/workspace/intergenic_up_downstream_faststructure/result/
seq 2 10 | parallel -j 10 "structure.py -K {} \
    --input=/home/vensin/workspace/intergenic_up_downstream_faststructure/194samples_filtered.intergenic_up_downstream.LD.prune  \
    --output=/home/vensin/workspace/population_structure/faststructure/result/LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > LD_faststructure_K_{}.log 2>&1" &

chooseK.py --input LD_nomissing_faststructure_K

## Model complexity that maximizes marginal likelihood = 
## Model components used to explain structure in data = 
