cd /home/vensin/workspace/Annovar

##  awk '$1 == "intergenic" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic.snp.id
##  6075980 intergenic.snp.id

awk '$1 == "intergenic" || $1 == "upstream" || $1 == "downstream" || $1 == "upstream;downstream" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic_up_downstream.snp.id
##  wc -l intergenic_up_downstream.snp.id
##  8210383 intergenic_up_downstream.snp.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic_up_downstream.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic_up_downstream.snp
##  After filtering, kept 8073005 out of a possible 11452745 Sites
##  Run Time = 2179.00 seconds

plink --vcf 194samples_filtered.intergenic_up_downstream.snp.recode.vcf \
    --indep-pairwise 50 5 0.2 \
    --out intergenic_up_downstream.LD  \
    --allow-extra-chr  \
    --set-missing-var-ids @:# 
## Pruning complete.  6958009 of 8073005 variants removed.
## 1114996

sed 's/:/ /g' intergenic_up_downstream.LD.prune.in > intergenic_up_downstream.LD.prune.in.snp.id

vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic_up_downstream.LD.prune.in.snp.id \
    --recode --recode-INFO-all \
    --out 194samples_filtered.intergenic_up_downstream.LD.prune
##  After filtering, kept 1114996 out of a possible 11452745 Sites
##  Run Time = 335.00 seconds

mkdir -p /home/vensin/workspace/intergenic_up_downstream_faststructure/result
cd /home/vensin/workspace/intergenic_up_downstream_faststructure

## vcf转bed格式
plink --vcf /home/vensin/workspace/Annovar/194samples_filtered.intergenic_up_downstream.LD.prune.recode.vcf \
    --make-bed   --keep-allele-order  --allow-extra-chr \
    --out /home/vensin/workspace/intergenic_up_downstream_faststructure/194samples_filtered.intergenic_up_downstream.LD.prune 

cd /home/vensin/workspace/intergenic_up_downstream_faststructure/result/
seq 2 10 | parallel -j 10 "structure.py -K {} \
    --input=/home/vensin/workspace/intergenic_up_downstream_faststructure/194samples_filtered.intergenic_up_downstream.LD.prune  \
    --output=/home/vensin/workspace/intergenic_up_downstream_faststructure/result/intergenic_up_downstream_LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > intergenic_up_downstream_LD_faststructure_K_{}.log 2>&1" &

chooseK.py --input intergenic_up_downstream_LD_faststructure_K

## Model complexity that maximizes marginal likelihood = 3
## Model components used to explain structure in data = 5

################################## /home/vensin/workspace/intergenic_up_downstream_faststructure/result




