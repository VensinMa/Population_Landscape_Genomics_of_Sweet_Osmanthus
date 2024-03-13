# python版本
## (bio-env) root@DESKTOP-O94TLKA:~/workspace/pop/vcf2phylip-2.8# python --version
## Python 3.10.13

# vcf2phylip.py版本
## (bio-env) root@DESKTOP-O94TLKA:~/workspace/pop# python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py  -v
## vcf2phylip.py 2.8

# treebest版本
## (bio-env) root@DESKTOP-O94TLKA:~/workspace/pop# treebest
## Program: TreeBeST (gene Tree Building guided by Species Tree)
## Version: 1.9.2 build 15May2023
## Contact: Heng Li <lh3@sanger.ac.uk>


###  以下为本地wsl运行

cd /root/workspace/pop
wget https://github.com/edgardomortiz/vcf2phylip/archive/refs/tags/v2.8.tar.gz
tar -zxvf v2.8.tar.gz


python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py --input 224_filtered.LD.pruned.noContig.recode.vcf --fasta --output-prefix 224_filtered.LD.pruned.noContig

treebest nj -b 1000 224_filtered.LD.pruned.noContig.min4.fasta > 224_filtered.LD.pruned.noContig.treebest.out
treebest nj -W -b 1000 224.filtered.LD.pruned.noContig.min4.fasta > 224.filtered.LD.pruned.noContig.treebest.unroot.tree.out &

sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 224_filtered.LD.pruned.noContig.treebest.out | awk '{printf $0}' > 224_filtered.LD.pruned.noContig.treebest.nwk
