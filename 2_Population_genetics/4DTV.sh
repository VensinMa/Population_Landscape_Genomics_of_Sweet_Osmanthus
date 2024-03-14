## (base) [guop@node01 229sample_LD_4DTV]$ cd /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV
## (base) [guop@node01 229sample_LD_4DTV]$ pwd
## /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV
## (base) [guop@node01 229sample_LD_4DTV]$  grep -cv "^#" 229.snpEffAnno.4dtv.vcf
## 18987

####################################################   PCA  ##################################################
## vcf转bed
plink --vcf 229.snpEffAnno.4dtv.vcf  --make-bed   --out 229_filtered.LD.4dtv.noContig  --keep-allele-order  --allow-extra-chr 

## 生成grm文件
gcta64 --bfile 229_filtered.LD.4dtv.noContig  --autosome  --make-grm  --out GA

## 进行PCA分析
gcta64 --grm GA --pca 130  --out 229_filtered.LD.4dtv.noContig_PCA_out

###################################################  NJ tree ##################################################

python  /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py  --input 229.snpEffAnno.4dtv.vcf  --fasta --output-prefix  229_filtered.LD.4dtv.noContig

treebest nj -b 1000 229_filtered.LD.4dtv.noContig.min4.fasta > 229_filtered.LD.4dtv.noContig.treebest.out

sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 229_filtered.LD.4dtv.noContig.treebest.out | awk '{printf $0}' > 229_filtered.LD.4dtv.noContig.treebest.nwk

###################################################  ML tree ##################################################

iqtree -s 229_filtered.LD.4dtv.noContig.min4.fasta  -m MFP  -B 1000  -keep-ident  -T 20 &









