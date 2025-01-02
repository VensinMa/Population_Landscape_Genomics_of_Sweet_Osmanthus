### 以下分析在 ubuntu24 系统进行

##  annovar官方网站的的下载需要使用邮箱注册后才可下载
# https://www.openbioinformatics.org/annovar/annovar_download_form.php

################################# annovar 安装 ##############################################
mkdir -p /home/vensin/workspace/Annovar/result && cd /home/vensin/workspace/Annovar
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

  ### NOTICE: Finished reading 11452900 lines from VCF file
  ### NOTICE: A total of 11452745 locus in VCF file passed QC threshold, representing 11452745 SNPs (8458075 transitions and 2994670 transversions) and 0 indels/substitutions
  ### NOTICE: Finished writing allele frequencies based on 2221832530 SNP genotypes (1640866550 transitions and 580965980 transversions) and 0 indels/substitutions for 194 samples

## 进行变异注释 (如果需要所有信息，添加-separate) -separate 将每种类型的变异分开注释到不同的文件中
annotate_variation.pl -geneanno --neargene 2000 -buildver  LYG.hic -dbtype refGene -outfile LYG.hic.snp.annovar -exonsort 194samples_filtered.annovar.input  ./
  ### NOTICE: Output files are written to LYG.hic.snp.annovar.variant_function, LYG.hic.snp.annovar.exonic_variant_function
  ### NOTICE: Reading gene annotation from ./LYG.hic_refGene.txt ... Done with 56392 transcripts (including 0 without coding sequence annotation) for 41252 unique genes
  ### NOTICE: Processing next batch with 5000000 unique variants in 5000000 input lines
  ### NOTICE: Finished analyzing 1000000 query variants
  ### NOTICE: Finished analyzing 2000000 query variants
  ### NOTICE: Finished analyzing 3000000 query variants
  ### NOTICE: Finished analyzing 4000000 query variants
  ### NOTICE: Reading FASTA sequences from ./LYG.hic_refGeneMrna.fa ... Done with 22635 sequences
  ### NOTICE: Processing next batch with 5000000 unique variants in 5000000 input lines
  ### NOTICE: Finished analyzing 1000000 query variants
  ### NOTICE: Finished analyzing 2000000 query variants
  ### NOTICE: Finished analyzing 3000000 query variants
  ### NOTICE: Finished analyzing 4000000 query variants
  ### NOTICE: Reading FASTA sequences from ./LYG.hic_refGeneMrna.fa ... Done with 21094 sequences
  ### NOTICE: Processing next batch with 1452745 unique variants in 1452745 input lines
  ### NOTICE: Finished analyzing 1000000 query variants
  ### NOTICE: Reading FASTA sequences from ./LYG.hic_refGeneMrna.fa ... Done with 5413 sequences
  


####################################   统计每种类型（突变位置）SNP的数量  ###############################################
# cat LYG.hic.snp.annovar.variant_function | cut -f 1 | sed 's/;/\n/g' | sort | uniq -c
cat LYG.hic.snp.annovar.variant_function | cut -f 1 | sort | uniq -c

 910824 downstream
 533562 exonic
    156 exonic;splicing
6075980 intergenic
2529888 intronic
   3436 splicing
1086201 upstream
 137378 upstream;downstream
 100575 UTR3
  74075 UTR5
    670 UTR5;UTR3

    781 downstream
    970 exonic
   3793 intergenic
   3120 intronic
      2 splicing
    625 upstream
     72 upstream;downstream
    193 UTR3
    118 UTR5
#######################################   统计外显子区域不同突变类型SNP的数量  #############################################
cat LYG.hic.snp.annovar.exonic_variant_function | awk '{print $2}' | sort | uniq -c

 290892 nonsynonymous
   6097 stopgain
    806 stoploss
 235923 synonymous


    511 nonsynonymous
      4 stopgain
    455 synonymous

 
