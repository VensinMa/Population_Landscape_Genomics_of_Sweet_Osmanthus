### 学院服务器上
#  /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.gatk.vcf

cd  /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools

###  剔除外类群样本   O_DSMX  O_MXL   O_XYWJM   保留浙南木犀 网脉木犀 蒙自桂花 O_ZNMX  NFM_1   O_MZGH
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.gatk.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_3_outgroup.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/197samples_filtered_renamed.snp.unanchor.final.gatk.vcf &

## After filtering, kept 197 out of 200 Individuals
## Outputting VCF file...
## After filtering, kept 141490137 out of a possible 141490137 Sites
## Run Time = 384836.00 seconds

vcftools --vcf 197samples_filtered_renamed.snp.unanchor.final.gatk.vcf.recode.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --minGQ 10 \
    --minQ 30 \
    --min-meanDP 6 \
    --max-missing 1 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out 197samples_filtered.snp.nomissing

###  剔除外类群样本   O_DSMX  O_MXL  O_XYWJM  O_MZGH  保留浙南木犀 网脉木犀  O_ZNMX  NFM_1  
vcftools --vcf /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/200samples_filtered_renamed.snp.unanchor.final.gatk.vcf \
    --recode --recode-INFO-all  --remove /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/removed_4_outgroup.id \
    --out /public1/guop/mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/vcftools/196samples_filtered_renamed.snp.unanchor.final.gatk.vcf &

## After filtering, kept 196 out of 200 Individuals
## Outputting VCF file...
## After filtering, kept 141490137 out of a possible 141490137 Sites
## Run Time = 384708.00 seconds

vcftools --vcf 196samples_filtered_renamed.snp.unanchor.final.gatk.vcf.recode.vcf \
    --min-alleles 2 \
    --max-alleles 2 \
    --minDP 5 \
    --minGQ 10 \
    --minQ 30 \
    --min-meanDP 6 \
    --max-missing 1 \
    --maf 0.05 \
    --recode \
    --recode-INFO-all \
    --out 196samples_filtered.snp.nomissing

# 1、将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/est-sfs/prepare_est-sfs
python vcf_to_estsfs.py  196samples_filtered.snp.nomissing.recode.vcf  DRS-7 LCJ-7
## 194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.recode_estsfs_input.txt

# 2、运行 est-sfs
cd /home/vensin/workspace/est-sfs/
est-sfs config-2outgroup.txt  /home/vensin/workspace/est-sfs/prepare_est-sfs/194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.plink_estsfs_input.txt  seedfile.txt output_file_sfs.txt  output_file_p_anc.txt

# 3、极性化原vcf文件
cd /home/vensin/workspace/est-sfs/
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.plink.vcf  output_file_p_anc.txt

# 4、计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
## pip3 install vcfpy  # 使用脚本前需要import vcf

