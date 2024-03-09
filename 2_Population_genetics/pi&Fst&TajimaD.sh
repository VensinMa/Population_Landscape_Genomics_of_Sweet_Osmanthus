
# 计算π值
vcftools --vcf input.vcf --window-pi 100000 --out pi_100kb

# 计算Tajima's D值
vcftools --vcf input.vcf --TajimaD 100000 --out tajimasD_100kb

# 计算FST
vcftools --vcf input.vcf --weir-fst-pop pop1.txt --weir-fst-pop pop2.txt --fst-window-size 100000 --fst-window-step 1000 --out fst_100kb_1kb
