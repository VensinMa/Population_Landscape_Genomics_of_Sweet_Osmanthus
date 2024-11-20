####  https://blog.csdn.net/2302_79242191/article/details/134630776?utm_medium=distribute.pc_relevant.none-task-blog-2~default~baidujs_baidulandingword~default-0-134630776-blog-105530364.235^v43^pc_blog_bottom_relevance_base2&spm=1001.2101.3001.4242.1&utm_relevant_index=1

### 0.准备
cd /home/vensin/software/snpEff_latest_core/snpEff  && chmod +x SnpSift.jar  && chmod +x snpEff.jar
echo "guihua.genome:guihua" >> /home/vensin/software/snpEff_latest_core/snpEff/snpEff.conf

##  mkdir -p data/XXX  
mkdir -p /home/vensin/software/snpEff_latest_core/snpEff/data/guihua 

## 在data下面存放两个文件  LYG.hic.fasta: 参考基因组  LYG.hic.gff: 注释文件，GFF3格式（也可以是GFF2格式）
cd /home/vensin/software/snpEff_latest_core/snpEff/data/

## java -jar snpEff.jar build -gff3 -v XXX

### 1.建库
java -jar /home/vensin/software/snpEff_latest_core/snpEff/snpEff.jar build -gff3 -v guihua -noCheckProtein -noCheckCds


### 2.注释  
java -jar /home/vensin/software/snpEff_latest_core/snpEff/snpEff.jar -c  /home/vensin/software/snpEff_latest_core/snpEff/snpEff.config \
  -ud 2000 -csvStats guihua.csv -htmlStats guihua.html -o vcf  guihua \
  /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink.polarized.vcf > /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff.vcf

java -jar /home/vensin/software/snpEff_latest_core/snpEff/snpEff.jar -c  /home/vensin/software/snpEff_latest_core/snpEff/snpEff.config \
  -ud 2000 -csvStats guihua.lof.csv -htmlStats guihua.lof.html -o vcf  -lof guihua \
  /home/vensin/workspace/est-sfs/prepare_est-sfs/197samples_filtered_3_outgroup.snp.nomissing.rename.plink.polarized.vcf > /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff.vcf

### 
grep '^#' /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff.vcf  > /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff_LOF.vcf && \
  grep 'LOF' /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff.vcf >> /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff_LOF.vcf

python /home/vensin/workspace/est-sfs/prepare_est-sfs/indv_GT_stats.py /home/vensin/workspace/est-sfs/197samples_filtered_3_outgroup.polarized.snpeff_LOF.vcf --output /home/vensin/workspace/est-sfs/indv_GT_stats_res.197samples_filtered_3_outgroup.polarized.snpeff_LOF.txt

python /home/vensin/workspace/est-sfs/prepare_est-sfs/indv_GT_stats.py /home/vensin/workspace/sift-lyjg/deleterious_197.vcf --output /home/vensin/workspace/sift-lyjg/indv_GT_stats_res.deleterious_197.txt
