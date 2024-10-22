##  annovar官方网站的的下载需要使用邮箱注册后才可下载
# https://www.openbioinformatics.org/annovar/annovar_download_form.php

################################# annovar 安装 ##############################################
mkdir -p /home/vensin/workspace/Annovar && cd /home/vensin/workspace/Annovar
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz   # wget https://mawenxin.cn/annovar.latest.tar.gz

tar -zxvf annovar.latest.tar.gz

## 安装注释文件格式转换工具 gffread gff3ToGenePred gtfToGenePred
conda install -c bioconda gffread

conda create -n ucsc
conda activate ucsc
conda install bioconda::ucsc-gff3togenepred
conda install bioconda::ucsc-gtftogenepred

## 添加环境变量 export PATH
echo 'export PATH=/home/vensin/workspace/Annovar/annovar:$PATH' >> ~/.bashrc
echo 'export PATH=/home/vensin/anaconda3/envs/ucsc/bin/:$PATH' >> ~/.bashrc
source ~/.bashrc

########################  gff格式转gtf格式  ############################################################################
gffread  /home/vensin/workspace/Annovar/final.gff -T -o LYG.hic.gtf

gtfToGenePred  -genePredExt LYG.hic.gtf LYG.hic_refGene.txt

# gff3ToGenePred  Of.genome.gff3  Of.genome_refGene.txt
# gff3文件开头必须是##gff-version 3

retrieve_seq_from_fasta.pl --format refGene --seqfile /home/vensin/workspace/Annovar/LYG.hic.fasta  LYG.hic_refGene.txt --out LYG.hic_refGeneMrna.fa

#############################  生成表格格式输入文件  ###################################################################
convert2annovar.pl -format vcf4 -allsample -withfreq \
  /home/vensin/workspace/Annovar/194samples_snp.nounanchor.renamed.filtered.vcftools.recode.vcf \
  > 194samples_filtered.annovar.input

## 进行变异注释 (如果需要所有信息，添加-separate) -separate 将每种类型的变异分开注释到不同的文件中
annotate_variation.pl -geneanno --neargene 2000 -buildver  LYG.hic -dbtype refGene -outfile LYG.hic.snp.annovar -exonsort 194samples_filtered.annovar.input  ./

##  NOTICE: Finished reading 12551354 lines from VCF file
##  NOTICE: A total of 12551267 locus in VCF file passed QC threshold, representing 12551267 SNPs (9218405 transitions and 3332862 transversions) and 0 indels/substitutions
##  NOTICE: Finished writing allele frequencies based on 2334535662 SNP genotypes (1714623330 transitions and 619912332 transversions) and 0 indels/substitutions for 186 samples

####################################   统计每种类型（突变位置）SNP的数量  ###############################################
cat AllSNP.annovar.variant_function | cut -f 1 | sed 's/;/\n/g' | sort | uniq -c
1253953 downstream
 592470 exonic
7548615 intergenic
1906150 intronic
   3022 splicing
1398247 upstream

#######################################   统计外显子区域不同突变类型SNP的数量  #############################################
cat AllSNP.annovar.exonic_variant_function | awk '{print $2}' | sort | uniq -c
 322622 nonsynonymous
   6955 stopgain
    769 stoploss
 261796 synonymous
    328 unknown
