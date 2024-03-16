cd /root/workspace/186sample
mkdir stat

# Getting allele frequency  使用 --freq 计算每个位点的等位基因频率
vcftools --vcf 224_filtered_vcftools.vcf.gz.recode.vcf  --freq2 --out 224_filtered_site_allele_frequency
vcftools --vcf 224_filtered.LD.pruned.recode.vcf  --freq2 --out 224_filtered.LD.pruned_site_allele_frequency

# Getting sequencing depth information 使用 --depth 计算每个个体的平均测序深度
vcftools --vcf  224_filtered_vcftools.vcf.gz.recode.vcf --depth --out 224_filtered_ind_mean_depth
vcftools --vcf  224_filtered.LD.pruned.recode.vcf --depth --out 224_filtered.LD.pruned_ind_mean_depth

# 使用 --site-mean-depth 计算每个变异位点的平均深度
vcftools --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --site-mean-depth --out 224_filtered_site_mean_depth
vcftools --vcf 224_filtered.LD.pruned.recode.vcf --site-mean-depth --out 224_filtered.LD.pruned_site_mean_depth

# 计算位点平均质量值
vcftools --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --site-quality --out 224_filtered_mean_quality
vcftools --vcf 224_filtered.LD.pruned.recode.vcf --site-quality --out 224_filtered.LD.pruned_site_mean_quality

# 计算每个个体数据缺失率
vcftools --vcf 224_filtered_vcftools.vcf.gz.recode.vcf --missing-indv --out 224_filtered_ind_missingindv
vcftools --vcf 224_filtered.LD.pruned.recode.vcf --missing-indv --out 224_filtered.LD.pruned_ind_missingindv
