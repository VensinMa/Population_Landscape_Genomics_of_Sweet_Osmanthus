cd  /public1/guop/mawx/workspace/186sample/fastastructure
mkdir fastastructure2_20_result
seq 2 20 | parallel -j 20 "structure.py -K {} --input=/public1/guop/mawx/workspace/186sample/186_filtered.LD.pruned.noContig --output=/public1/guop/mawx/workspace/186sample/fastastructure/fastastructure2_20_result/LD_faststructure_K_{} --cv=5 --prior=simple --seed=123 > LD_faststructure_K_{}.log 2>&1" &
