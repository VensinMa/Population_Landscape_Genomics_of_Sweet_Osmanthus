##  含外类群构树
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/ && mkdir NJtree

###  剔除外类群样本 NFM_1   O_DSMX  O_MXL   O_MZGH  O_XYWJM  只保留浙南木犀 O_ZNMX
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_5_outgroup_keep_ZNMX.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/195samples_filtered.LD.pruned.removed_5_outgroup_keep_ZNMX &

###  剔除外类群样本 NFM_1   O_MXL   O_MZGH  O_XYWJM  O_ZNMX  只保留短丝木犀 O_DSMX
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_5_outgroup_keep_DSMX.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/195samples_filtered.LD.pruned.removed_5_outgroup_keep_DSMX &

###  剔除外类群样本 NFM_1  O_DSMX  O_MXL   O_XYWJM  O_ZNMX   只保留蒙自桂花 O_MZGH
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_5_outgroup_keep_MZGH.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/195samples_filtered.LD.pruned.removed_5_outgroup_keep_MZGH &

###  剔除外类群样本  O_DSMX  O_MXL   O_MZGH  O_XYWJM  O_ZNMX 只保留网脉木犀 NFM_1
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_5_outgroup_keep_WMMX.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/195samples_filtered.LD.pruned.removed_5_outgroup_keep_WMMX &


    
#######################      VCF2Dis 计算距离矩阵构树   ##################################
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree
VCF2Dis -i  195samples_filtered.LD.pruned.removed_5_outgroup_keep_MZGH.recode.vcf  -o 195samples_filtered.LD.pruned.removed_5_outgroup_keep_MZGH_dis.mat    -Rand 0.25 &
VCF2Dis -i  195samples_filtered.LD.pruned.removed_5_outgroup_keep_DSMX.recode.vcf  -o 195samples_filtered.LD.pruned.removed_5_outgroup_keep_DSMX_dis.mat    -Rand 0.25 &
VCF2Dis -i  195samples_filtered.LD.pruned.removed_5_outgroup_keep_WMMX.recode.vcf  -o 195samples_filtered.LD.pruned.removed_5_outgroup_keep_WMMX_dis.mat    -Rand 0.25 &
VCF2Dis -i  195samples_filtered.LD.pruned.removed_5_outgroup_keep_ZNMX.recode.vcf  -o 195samples_filtered.LD.pruned.removed_5_outgroup_keep_ZNMX_dis.mat    -Rand 0.25 &

VCF2Dis -i  /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    -o /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/200samples_filtered.LD.pruned_dis.mat  \
    -Rand 0.25 &


#######################     treebest  构树   ###################################
cd /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree
python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 195samples_filtered.LD.pruned.removed_5_outgroup_keep_MZGH.recode.vcf --fasta --output-prefix 195samples_filtered.LD.pruned.removed_5_outgroup_keep_MZGH
python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 195samples_filtered.LD.pruned.removed_5_outgroup_keep_DSMX.recode.vcf --fasta --output-prefix 195samples_filtered.LD.pruned.removed_5_outgroup_keep_DSMX
python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 195samples_filtered.LD.pruned.removed_5_outgroup_keep_WMMX.recode.vcf --fasta --output-prefix 195samples_filtered.LD.pruned.removed_5_outgroup_keep_WMMX
python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 195samples_filtered.LD.pruned.removed_5_outgroup_keep_ZNMX.recode.vcf --fasta --output-prefix 195samples_filtered.LD.pruned.removed_5_outgroup_keep_ZNMX


python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf \
    --fasta --output-prefix /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/NJtree/200samples_filtered.LD.pruned
treebest nj -b 1000 200samples_filtered.LD.pruned.min4.fasta > 200samples_filtered.LD.pruned.treebest.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 200samples_filtered.LD.pruned.treebest.out | awk '{printf $0}' > 200samples_filtered.LD.pruned.treebest.nwk










