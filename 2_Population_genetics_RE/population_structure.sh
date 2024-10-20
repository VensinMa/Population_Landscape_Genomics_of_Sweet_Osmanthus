###  1、fastastructure

mkdir -p /home/vensin/workspace/population_structure/fastastructure/result
cd /home/vensin/workspace/population_structure
## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
    --make-bed   --out /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned  --keep-allele-order  --allow-extra-chr  

cd /home/vensin/workspace/population_structure/fastastructure/result/
seq 2 20 | parallel -j 8 "structure.py -K {} \
    --input=/home/vensin/workspace/population_structure/194samples_filtered.LD.pruned  \
    --output=/home/vensin/workspace/population_structure/fastastructure/result/LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > LD_faststructure_K_{}.log 2>&1" &


###  2、admixture
mkdir -p /home/vensin/workspace/population_structure/admixture/result
##  更改染色体名称
cd /home/vensin/workspace/population_structure/admixture
sed -E 's/Superscaffold([0-9]{1})\b/Chr0\1/; s/Superscaffold([1-9][0-9]*)/Chr\1/' \
    /home/vensin/workspace/population_structure/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
    > /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.recode.vcf
    
## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.recode.vcf \
    --make-bed   --out /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr  --keep-allele-order  --allow-extra-chr  
    
seq 2 20 | parallel -j 8 'admixture /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.bed {} --cv | tee ./result/admixture_K{}.log' &

seq 20 -1 8 | parallel -j 8 'admixture /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.bed {} --cv | tee ./result/admixture_K{}.log' &




