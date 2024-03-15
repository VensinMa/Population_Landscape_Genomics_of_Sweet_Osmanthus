## (base) [guop@node01 229sample_LD_4DTV]$ cd /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV
## (base) [guop@node01 229sample_LD_4DTV]$ pwd
## /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV
## (base) [guop@node01 229sample_LD_4DTV]$  grep -cv "^#" 229.snpEffAnno.4dtv.vcf
## 18987

####################################################   PCA  ##################################################
## vcf转bed
plink --vcf /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV/229.snpEffAnno.4dtv.vcf  --make-bed   --out 229_filtered.LD.4dtv.noContig  --keep-allele-order  --allow-extra-chr 

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

#################################################  faststructure  ##############################################
cd /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV/structure/faststructure
plink --vcf /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV/229.snpEffAnno.4dtv.vcf  --make-bed   --out 229_filtered.LD.4dtv.noContig  --keep-allele-order  --allow-extra-chr  
mkdir fast_result
seq 2 30 | parallel -j 30 "structure.py -K {} --input=229_filtered.LD.4dtv.noContig --output=fast_result/229_filtered.LD.4dtv.noContig_faststructure_{} --cv=5 --prior=logistic --seed=123 > faststructure.{}.log 2>&1" &

#################################################  admixture  ##############################################
cd /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV/structure/admixture
## 需要bed格式文件 
plink --vcf /public1/guop/mawx/workspace/rxg_wild_snpcalling/pop/229sample_LD_4DTV/229.snpEffAnno.4dtv.vcf  --make-bed   --out 229_filtered.LD.4dtv.noContig  --keep-allele-order  --allow-extra-chr  

# 计算K=2到30
seq 2 30 | parallel -j 30 "admixture --cv  229_filtered.LD.4dtv.noContig.bed {} 1>admix.{}.log 2>&1" &

# 确定具有最小CV值的K值为最佳分群数
cat *.log | grep "CV"

mkdir result
cp  ./*.Q result/









