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
python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py --input 229_filtered.LD.pruned.noContig.recode.vcf --fasta --output-prefix 229_filtered.LD.pruned.noContig
python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py --input modi_224.snpEffAnno.4dtv.nocontig.vcf   --fasta --output-prefix 224.snpEffAnno.4dtv.nocontig

/root/workspace/pop/VCF2Dis-1.50/bin/VCF2Dis -i  229_filtered.LD.pruned.noContig.recode.vcf -o 229_filtered.LD.pruned.noContig_dis.mat

treebest nj -b 1000 224_filtered.LD.pruned.noContig.min4.fasta > 224_filtered.LD.pruned.noContig.treebest.out
treebest nj -W -b 1000 224.filtered.LD.pruned.noContig.min4.fasta > 224.filtered.LD.pruned.noContig.treebest.unroot.tree.out &
treebest nj -W -b 1000 224.snpEffAnno.4dtv.nocontig.min4.fasta > 224.snpEffAnno.4dtv.nocontig.treebest.unroot.tree.out 

sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 224_filtered.LD.pruned.noContig.treebest.out | awk '{printf $0}' > 224_filtered.LD.pruned.noContig.treebest.nwk
224.snpEffAnno.4dtv.nocontig.treebest.unroot.tree.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 224.snpEffAnno.4dtv.nocontig.treebest.unroot.tree.out | awk '{printf $0}' > 224.snpEffAnno.4dtv.nocontig.treebest.unroot.tree.nwk





