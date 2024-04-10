cd ~/workspace/186sample/stat

# 计算π值
vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --window-pi 100000  --window-pi-step 1000 --out LD.pruned.noContig_pi_100kb
##  After filtering, kept 186 out of 186 Individuals
##  Outputting Windowed Nucleotide Diversity Statistics...
##  After filtering, kept 1789805 out of a possible 1789805 Sites
##  Run Time = 20.00 seconds
vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --window-pi 100000 --out LD.pruned.noContig_pi_100kb_nostep

# 计算Tajima's D值
vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --TajimaD 100000 --out LD.pruned.noContig_tajimasD_100kb
##  After filtering, kept 186 out of 186 Individuals
##  Outputting Tajima's D Statistic...
##  After filtering, kept 1789805 out of a possible 1789805 Sites
##  Run Time = 18.00 seconds

vcftools --vcf ../186_filtered.LD.pruned.noContig.recode.vcf --window-pi 100000  --window-pi-step 1000 --out LD.pruned.noContig_pi_100kb
awk '{$4=""; print $0}' LD.pruned.noContig_pi_100kb_nostep.windowed.pi >  LD.pruned.noContig_pi_100kb_nostep.windowed.pi.txt
awk '{$4=""; print $0}' LD.pruned.noContig_pi_100kb_nostep.windowed.pi >  LD.pruned.noContig_pi_100kb_nostep.windowed.pi.txt



###############################   LD位点 计算pi   #######################################################
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central-East_samples.txt --window-pi 100000 --out Central-East_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep West_samples.txt --window-pi 100000 --out West_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central_samples.txt --window-pi 100000 --out Central_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep East_samples.txt --window-pi 100000 --out East_samples_pi


(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central-East_samples_pi.windowed.pi
0.000527711
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' West_samples_pi.windowed.pi
0.000353991
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central_samples_pi.windowed.pi
0.000484129
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' East_samples_pi.windowed.pi
0.000529708



##############################  all 高质量SNP 计算pi ###########################################################
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central-East_samples_allsnp_pi.windowed.pi
0.00388759
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' West_samples_allsnp_pi.windowed.pi
0.00268107



vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central-East_samples.txt --window-pi 100000 --out Central-East_samples_allsnp_pi
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep West_samples.txt --window-pi 100000 --out West_samples_allsnp_pi
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central_samples.txt --window-pi 100000 --out Central_samples_allsnp_pi 
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep East_samples.txt --window-pi 100000 --out East_samples_allsnp_pi

















