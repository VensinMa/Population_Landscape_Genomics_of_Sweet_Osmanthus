##  fastStructure使用的输入文件基本与admixture 一致，bed格式

cd /public1/guop/mawx/workspace/pop_genetic/faststructure
mkdir fast_result
seq 2 30 | parallel -j 30 "structure.py -K {} --input=/public1/guop/mawx/workspace/pop_genetic/224.filtered.LD.pruned.noContig --output=fast_result/224.filtered.LD.pruned.noContig_{} --cv=5 --prior=logistic --seed=100 > log.{}.txt 2>&1" &

