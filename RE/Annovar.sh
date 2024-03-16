##  annovar官方网站的的下载需要使用邮箱注册后才可下载
# http://download.openbioinformatics.org/annovar_download_form.php

################################# annovar 安装 ##############################################
cd /root/workspace/annovar
wget http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/annovar.latest.tar.gz
tar -zxvf annovar.latest.tar.gz

## 安装注释文件格式转换工具 gffread gff3ToGenePred gtfToGenePred
conda install -c bioconda gffread

conda create -n ucsc
conda activate ucsc
conda install bioconda::ucsc-gff3togenepred
conda install bioconda::ucsc-gtftogenepred

## 添加环境变量 export PATH
vi ~/.bashrc	
export PATH=/root/workspace/annovar/annovar:$PATH 
export PATH=~/anaconda2/envs/ucsc/bin/:$PATH
source ~/.bashrc

########################  gff格式转gtf格式  ######################################################
gffread  /root/workspace/genome/guihua.genomic.gff3 -T -o  guihua.genomic.gtf

gtfToGenePred  -genePredExt guihua.genomic.gtf guihua.genomic_refGene.txt

# gff3ToGenePred  Of.genome.gff3  Of.genome_refGene.txt
# gff3文件开头必须是##gff-version 3

retrieve_seq_from_fasta.pl --format refGene --seqfile /root/workspace/genome/Osmanthus.genomic.fasta  guihua.genomic_refGene.txt --out guihua.genomic_refGeneMrna.fa

#############################  生成表格格式输入文件  #################################################
convert2annovar.pl -format vcf4 -allsample -withfreq  /root/workspace/186sample/186_filtered_vcftools.noContig.recode.vcf > 186_filtered_vcftools.noContig.annovar.input

## 进行变异注释 (如果需要所有信息，添加-separate) -separate 将每种类型的变异分开注释到不同的文件中
annotate_variation.pl -geneanno --neargene 2000 -buildver  guihua.genomic -dbtype refGene -outfile AllSNP.annovar -exonsort 186_filtered_vcftools.noContig.annovar.input  ./










