cd /home/vensin/workspace/Annovar

awk '$1 == "intergenic" {print $3 " " $4}' LYG.hic.snp.annovar.variant_function > intergenic.snp.id
##  6075980 intergenic.snp.id

mkdir -p /home/vensin/workspace/intergenic_faststructure/result
vcftools --vcf 194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
    --positions intergenic.snp.id \
    --recode --recode-INFO-all \
    --out /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.snp
## After filtering, kept 6075980 out of a possible 11452745 Sites
## Run Time = 1632.00 seconds

plink --vcf /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.snp.recode.vcf \
    --indep-pairwise 50 5 0.2 \
    --out LD  \
    --allow-extra-chr  \
    --set-missing-var-ids @:# 
##  Pruning complete.  5219387 of 6075980 variants removed.
##  Marker lists written to LD.prune.in and LD.prune.out .

sed 's/:/ /g' LD.prune.in > intergenic.LD.prune.in.snp.id

vcftools --vcf /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.snp.recode.vcf \
    --positions /home/vensin/workspace/Annovar/intergenic.LD.prune.in.snp.id \
    --recode --recode-INFO-all \
    --out /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune
## After filtering, kept 856593 out of a possible 6075980 Sites
## Run Time = 246.00 seconds

cd /home/vensin/workspace/intergenic_faststructure

## vcf转bed格式
plink --vcf /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune.recode.vcf \
    --make-bed   --keep-allele-order  --allow-extra-chr \
    --out /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune 

cd /home/vensin/workspace/intergenic_faststructure/result/
seq 2 10 | parallel -j 10 "structure.py -K {} \
    --input=/home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune  \
    --output=/home/vensin/workspace/intergenic_faststructure/result/intergenic_LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > intergenic_LD_faststructure_K_{}.log 2>&1" &

chooseK.py --input intergenic_LD_faststructure_K

## Model complexity that maximizes marginal likelihood = 2
## Model components used to explain structure in data = 7

############################ admixture ##################################

mkdir /home/vensin/workspace/intergenic_admixture/result

##  更改染色体名称
cd /home/vensin/workspace/intergenic_admixture/
sed -E 's/Superscaffold([0-9]{1})\b/Chr0\1/; s/Superscaffold([1-9][0-9]*)/Chr\1/' \
    /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune.recode.vcf \
    > /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune.Superscaffold2Chr.recode.vcf
    
## vcf转bed格式
plink --vcf /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune.Superscaffold2Chr.recode.vcf \
    --make-bed   --out /home/vensin/workspace/intergenic_faststructure/194samples_filtered.intergenic.LD.prune.Superscaffold2Chr  --keep-allele-order  --allow-extra-chr  
    
seq 2 10 | parallel -j 10 'admixture /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.bed {} --cv | tee ./result/admixture_K{}.log' &
cd /home/vensin/workspace/intergenic_admixture/result
# 确定具有最小CV值的K值为最佳分群数
cat *.log | grep "CV"

