cd /root/workspace/186sample

plink --vcf 186_filtered_vcftools.noContig.recode.vcf --indep-pairwise 50 5 0.2 --out LD  --allow-extra-chr  --set-missing-var-ids @:#  
##  wc -l  LD.prune.in
##  LD.prune.in

sed 's/:/ /g' LD.prune.in > LD.prune.in.VCFTOOLS.txt

vcftools --vcf 186_filtered_vcftools.noContig.recode.vcf  --positions LD.prune.in.VCFTOOLS.txt --recode --recode-INFO-all --out 186_filtered.LD.pruned.noContig
##  After filtering, kept 186 out of 186 Individuals
##  Outputting VCF file...
##  After filtering, kept 1789805 out of a possible 12551267 Sites
##  Run Time = 331.00 seconds





plink --vcf 186_filtered.LD.pruned.noContig.recode.vcf --recode --out 186_filtered.LD.pruned.noContig.ped --allow-extra-chr  --set-missing-var-ids @:#  
plink --file 186_filtered.LD.pruned.noContig.ped --recodeA --out 186_filtered.LD.pruned.noContig.raw --allow-extra-chr  --set-missing-var-ids @:#  
