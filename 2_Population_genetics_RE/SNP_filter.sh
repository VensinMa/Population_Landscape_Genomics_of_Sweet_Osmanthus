## 剔除多余的网脉木犀 最后保留了6个outgroup个体 O-WMMX O-DSMX O-MXL O-MZGH O-XYWJM O-ZNMX 和 194个内类群个体
vcftools --vcf all.filtered.snp.unanchor.vcf --keep keep_200sample.id --recode --recode-INFO-all --out 200samples_filtered.snp.unanchor

bcftools reheader -s 200sample_rename.id  200samples_filtered.snp.unanchor.recode.vcf -o 200samples_filtered_renamed.snp.unanchor.recode.vcf 
