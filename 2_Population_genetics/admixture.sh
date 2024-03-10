# 需要bed格式文件
plink --vcf 224__filtered.LD.pruned.noContig.recode.vcf --make-bed   --out 224_filtered.LD.pruned.noContig --keep-allele-order --allow-extra-chr  

# 计算K=2到20
seq 2 20 | parallel -j 20 "admixture --cv  224_filtered.LD.pruned.noContig.bed {} 1>admix.{}.log 2>&1"

# 确定具有最小CV值的K值为最佳分群数
cat *.log | grep "CV"

mkdir result
cp  ./*.Q result/

Rscript  draw_admixture.R  result  224.filtered.LD.pruned.noContig.nosex   224.filtered.LD.pruned.noContig.admixture

Rscript.exe  draw_admixture.R result 224.filtered.LD.pruned.noContig.nosex  224.filtered.LD.pruned.noContig.admixture
