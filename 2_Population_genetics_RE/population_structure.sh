###  1、fastastructure             0.2 缺失

mkdir -p /home/vensin/workspace/population_structure/faststructure/result
cd /home/vensin/workspace/population_structure
## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
    --make-bed   --out /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned  --keep-allele-order  --allow-extra-chr  

cd /home/vensin/workspace/population_structure/faststructure/result/
seq 2 20 | parallel -j 8 "structure.py -K {} \
    --input=/home/vensin/workspace/population_structure/194samples_filtered.LD.pruned  \
    --output=/home/vensin/workspace/population_structure/faststructure/result/LD_faststructure_K \
    --cv=5 --prior=simple --seed=123 > LD_faststructure_K_{}.log 2>&1" &
cd /home/vensin/workspace/population_structure/faststructure/result
chooseK.py --input LD_nomissing_faststructure_K
## Model complexity that maximizes marginal likelihood = 4
## Model components used to explain structure in data = 7


###  2、admixture             0.2 缺失
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
cd /home/vensin/workspace/population_structure/admixture/result
# 确定具有最小CV值的K值为最佳分群数
cat *.log | grep "CV"
'''
CV error (K=10): 0.47088
CV error (K=11): 0.47451
CV error (K=12): 0.47303
CV error (K=13): 0.49311
CV error (K=14): 0.50333
CV error (K=15): 0.50351
CV error (K=16): 0.52020
CV error (K=17): 0.52149
CV error (K=18): 0.53089
CV error (K=19): 0.53875
CV error (K=20): 0.54829
CV error (K=2): 0.44811
CV error (K=3): 0.44328
CV error (K=4): 0.44051
CV error (K=5): 0.44031
CV error (K=6): 0.44193
CV error (K=7): 0.44250
CV error (K=8): 0.44496
CV error (K=9): 0.46047
'''


###  3、fastastructure             0缺失

mkdir -p /home/vensin/workspace/population_structure/faststructure/nomissing_result
cd /home/vensin/workspace/population_structure

vcftools --vcf /home/vensin/workspace/population_structure/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf \
    --max-missing 1 \
    --recode \
    --recode-INFO-all \
    --out /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned.nomissing
    
## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned.nomissing.recode.vcf \
    --make-bed   --out /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned.nomissing  --keep-allele-order  --allow-extra-chr  

cd /home/vensin/workspace/population_structure/faststructure/nomissing_result
seq 2 10 | parallel -j 10 "structure.py -K {} \
    --input=/home/vensin/workspace/population_structure/194samples_filtered.LD.pruned.nomissing  \
    --output=/home/vensin/workspace/population_structure/faststructure/result/LD_nomissing_faststructure_K \
    --cv=5 --prior=simple --seed=1234 > LD_nomissing_faststructure_K_{}.log 2>&1" &


###  4、admixture             0缺失
mkdir -p /home/vensin/workspace/population_structure/admixture/nomissing_result
##  更改染色体名称
cd /home/vensin/workspace/population_structure/admixture
sed -E 's/Superscaffold([0-9]{1})\b/Chr0\1/; s/Superscaffold([1-9][0-9]*)/Chr\1/' \
    /home/vensin/workspace/population_structure/194samples_filtered.LD.pruned.nomissing.recode.vcf \
    > /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.nomissing.Superscaffold2Chr.recode.vcf
    
## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.nomissing.Superscaffold2Chr.recode.vcf \
    --make-bed   --out /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.nomissing.Superscaffold2Chr  --keep-allele-order  --allow-extra-chr  

cd /home/vensin/workspace/population_structure/admixture/nomissing_result
seq 2 10 | parallel -j 10 'admixture /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.nomissing.Superscaffold2Chr.bed {} --cv | tee ./admixture_nomissing_K{}.log' &

seq 20 -1 8 | parallel -j 8 'admixture /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.nomissing.Superscaffold2Chr.bed {} --cv | tee ./admixture_nomissing_K{}.log' &



