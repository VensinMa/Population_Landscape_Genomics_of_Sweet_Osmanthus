cd ~/workspace/186sample/stat

# 计算π值
vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --window-pi 100000  --window-pi-step 1000 --out LD.pruned.noContig_pi_100kb
##  After filtering, kept 186 out of 186 Individuals
##  Outputting Windowed Nucleotide Diversity Statistics...
##  After filtering, kept 1789805 out of a possible 1789805 Sites
##  Run Time = 20.00 seconds


# 计算Tajima's D值
vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --TajimaD 100000 --out LD.pruned.noContig_tajimasD_100kb
##  After filtering, kept 186 out of 186 Individuals
##  Outputting Tajima's D Statistic...
##  After filtering, kept 1789805 out of a possible 1789805 Sites
##  Run Time = 18.00 seconds



vcftools --vcf ../186.snpEffAnno.4dtv.LD.vcf --window-pi 100000 --out LD.pruned.4DTV.noContig_pi_100kb




vcftools --vcf ../186.snpEffAnno.4dtv.LD.vcf --TajimaD 100000 --out LD.pruned.noContig_tajimasD_100kb








