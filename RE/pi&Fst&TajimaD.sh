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

## 计算 FST 值
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --weir-fst-pop  West_samples.txt --weir-fst-pop Central-East_samples.txt  --fst-window-size 10000 --out West_Central-East_FST_10k

###############################   LD位点 计算pi   #######################################################
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central-East_samples.txt --window-pi 100000 --out Central-East_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep West_samples.txt --window-pi 100000 --out West_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central_samples.txt --window-pi 100000 --out Central_samples_pi
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep East_samples.txt --window-pi 100000 --out East_samples_pi

(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' West_samples_pi.windowed.pi
0.000353991
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central-East_samples_pi.windowed.pi
0.000527711
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central_samples_pi.windowed.pi
0.000484129
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' East_samples_pi.windowed.pi
0.000529708
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' LD.pruned.noContig_pi_100kb.windowed.pi
0.000519577
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' LD.pruned.noContig_pi_100kb_nostep.windowed.pi
0.000519501

##############################  all 高质量SNP 计算pi ###########################################################
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central-East_samples.txt --window-pi 100000 --out Central-East_samples_allsnp_pi
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep West_samples.txt --window-pi 100000 --out West_samples_allsnp_pi
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central_samples.txt --window-pi 100000 --out Central_samples_allsnp_pi 
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep East_samples.txt --window-pi 100000 --out East_samples_allsnp_pi
vcftools --vcf /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --window-pi 100000 --out 186_filtered_vcftools.noContig_pi_100kb_nostep_allsnp_pi

(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central-East_samples_allsnp_pi.windowed.pi
0.00388759
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' West_samples_allsnp_pi.windowed.pi
0.00268107
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' Central_samples_allsnp_pi.windowed.pi
0.00355378
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }'  East_samples_allsnp_pi.windowed.pi
0.00388855
(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $5; count++ } END { print total/count }' 186_filtered_vcftools.noContig_pi_100kb_nostep_allsnp_pi.windowed.pi
0.00392481


###############################   LD位点 计算 TajimaD   #######################################################
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central-East_samples.txt --TajimaD 100000 --out Central-East_samples_TajimaD &
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep West_samples.txt --TajimaD 100000 --out West_samples_TajimaD &
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep Central_samples.txt --TajimaD 100000 --out Central_samples_TajimaD &
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf --keep East_samples.txt --TajimaD 100000 --out East_samples_TajimaD &
vcftools --vcf /root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf  --TajimaD 100000 --out East_samples_TajimaD &




awk 'NR>1 { total += $4; count++ } END { print total/count }' Central-East_samples_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' West_samples_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' Central_samples_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' East_samples_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' LD.pruned.noContig_tajimasD_100kb.Tajima.D
1.81442
0.335453
0.5834
1.7114
1.90324


(bio-env) root@DESKTOP-O94TLKA:~/workspace/186sample/stat# awk 'NR>1 { total += $4; count++ } END { print total/count }'  LD.pruned.noContig_tajimasD_100kb.Tajima.D
1.90324

##############################  all 高质量SNP 计算TajimaD ###########################################################
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central-East_samples.txt --TajimaD 100000 --out Central-East_samples_allsnp_pTajimaD & 
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep West_samples.txt --TajimaD 100000 --out West_samples_allsnp_TajimaD &
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep Central_samples.txt --TajimaD 100000 --out Central_samples_allsnp_TajimaD  &
vcftools --vcf  /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --keep East_samples.txt --TajimaD 100000 --out East_samples_allsnp_TajimaD &
vcftools --vcf /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf --TajimaD 100000 --out 186_filtered_vcftools.noContig_TajimaD_100kb_nostep_allsnp_TajimaD &

awk 'NR>1 { total += $4; count++ } END { print total/count }' Central-East_samples_allsnp_pTajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' West_samples_allsnp_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' Central_samples_allsnp_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' East_samples_allsnp_TajimaD.Tajima.D
awk 'NR>1 { total += $4; count++ } END { print total/count }' 186_filtered_vcftools.noContig_TajimaD_100kb_nostep_allsnp_TajimaD.Tajima.D
1.96362
0.586069
0.68638
1.84944
2.13404




















