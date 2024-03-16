cd /mnt/e/mwx/workspace/pop/
grep -v "^Contig" 186_filtered_vcftools.recode.vcf

plink --vcf /mnt/e/mwx/workspace/pop/186_filtered_vcftools.recode.vcf --indep-pairwise 50 5 0.2 --out LD  --allow-extra-chr  --set-missing-var-ids @:#  
##  wc -l  LD.prune.in
##  
grep -v "Contig" LD.prune.in > LD.prune.noContig.in
## wc -l LD.prune.noContig.in
## 1727654 LD.prune.noContig.in
