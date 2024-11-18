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
    
## After filtering, kept 197 out of 197 Individuals
## Outputting VCF file...
## After filtering, kept 996932 out of a possible 141490137 Sites
## Run Time = 43840.00 seconds

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

## After filtering, kept 196 out of 196 Individuals
## Outputting VCF file...
## After filtering, kept 1104492 out of a possible 141490137 Sites
## Run Time = 44133.00 seconds

### 本地ubuntu
cd /home/vensin/workspace/est-sfs/prepare_est-sfs 
bcftools reheader -s 200sample_rename.id  197samples_filtered.snp.nomissing.recode.vcf  -o 197samples_filtered_3_outgroup.snp.nomissing.rename.recode.vcf
bcftools reheader -s 200sample_rename.id  196samples_filtered.snp.nomissing.recode.vcf  -o 196samples_filtered_2_outgroup.snp.nomissing.rename.recode.vcf
plink --vcf 197samples_filtered_3_outgroup.snp.nomissing.rename.recode.vcf --out 197samples_filtered_3_outgroup.snp.nomissing.rename.plink --recode vcf-iid  --allow-extra-chr  --keep-allele-order  --set-missing-var-ids @:# 
plink --vcf 196samples_filtered_2_outgroup.snp.nomissing.rename.recode.vcf --out 196samples_filtered_2_outgroup.snp.nomissing.rename.plink --recode vcf-iid  --allow-extra-chr  --keep-allele-order  --set-missing-var-ids @:# 

# 1、将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/est-sfs/prepare_est-sfs
python vcf_to_estsfs.py  196samples_filtered_2_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX
python vcf_to_estsfs.py  197samples_filtered_3_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX O-MZGH

'''
(base) vensin@ubuntu24-04:~/workspace/est-sfs/prepare_est-sfs$ python vcf_to_estsfs.py  196samples_filtered_2_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX
VCF file has been converted to est-sfs input file.
Total sites processed: 1104492
Successfully kept sites: 948992
(base) vensin@ubuntu24-04:~/workspace/est-sfs/prepare_est-sfs$ python vcf_to_estsfs.py  197samples_filtered_3_outgroup.snp.nomissing.rename.plink.vcf  O-WMMX O-ZNMX O-MZGH
VCF file has been converted to est-sfs input file.
Total sites processed: 996932
Successfully kept sites: 808574
'''
'''
-rw-rw-r-- 1 vensin vensin  15M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_processing.log
-rw-rw-r-- 1 vensin vensin  22M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
-rw-rw-r-- 1 vensin vensin  25M 11月 18 19:20 196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs_input.txt
drwxrwxr-x 2 vensin vensin 4.0K 11月 18 19:20 ./
-rw-rw-r-- 1 vensin vensin  18M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_processing.log
-rw-rw-r-- 1 vensin vensin  19M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
-rw-rw-r-- 1 vensin vensin  28M 11月 18 19:20 197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs_input.txt
'''

# 2、运行 est-sfs
cd /home/vensin/workspace/est-sfs/
'''
(base) vensin@ubuntu24-04:~/workspace/est-sfs$ cat config-2outgroup.txt config-3outgroup.txt
n_outgroup 2
model 2
nrandom 100
n_outgroup 3
model 2
nrandom 100
'''
est-sfs config-2outgroup.txt  /home/vensin/workspace/est-sfs/prepare_est-sfs/196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs_input.txt  seedfile.txt 2_outgroup_output_file_sfs.txt  2_outgroup_output_file_p_anc.txt
est-sfs config-3outgroup.txt  /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs_input.txt  seedfile.txt 3_outgroup_output_file_sfs.txt  3_outgroup_output_file_p_anc.txt

# 3、极性化原vcf文件
cd /home/vensin/workspace/est-sfs/
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/196samples_filtered_2_outgroup.snp.nomissing.rename.plink.vcf  2_outgroup_output_file_p_anc.txt /home/vensin/workspace/est-sfs/prepare_est-sfs/196samples_filtered_2_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
# Total polarized sites: 948992
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink.vcf  3_outgroup_output_file_p_anc.txt /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink_estsfs.positions.txt
# Total polarized sites: 808574

# 4、计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
python indv_GT_stats.py 197samples_filtered_3_outgroup.snp.nomissing.rename.plink.polarized.vcf 





