# 1、将vcf文件转换为est-sfs软件所需的输入文件格式
cd /home/vensin/workspace/est-sfs/prepare_est-sfs
python vcf_to_estsfs.py  194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.plink.vcf  DRS-7 LCJ-7
## 194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.recode_estsfs_input.txt

# 2、运行 est-sfs
cd /home/vensin/workspace/est-sfs/
est-sfs config-2outgroup.txt  /home/vensin/workspace/est-sfs/prepare_est-sfs/194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.plink_estsfs_input.txt  seedfile.txt output_file_sfs.txt  output_file_p_anc.txt

# 3、极性化原vcf文件
cd /home/vensin/workspace/est-sfs/
python vcf_polarize.py /home/vensin/workspace/est-sfs/prepare_est-sfs/194samples_filtered.intergenic_up_downstream.LD.prune.nomissing.plink.vcf  output_file_p_anc.txt

# 4、计算内类群个体突变的基因型数量 （相对于新的参考基因型 —— 祖先等位基因）
## pip3 install vcfpy  # 使用脚本前需要import vcf

