cd /root/workspace/186sample

python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py --input 186_filtered.LD.pruned.noContig.recode.vcf --fasta --output-prefix 186_filtered.LD.pruned.noContig
python /root/workspace/pop/vcf2phylip-2.8/vcf2phylip.py --input 186.snpEffAnno.4dtv.LD.vcf --fasta --output-prefix 186_filtered.LD.pruned.4DTV.noContig

treebest nj -b 1000 186_filtered.LD.pruned.4DTV.noContig.min4.fasta > 186_filtered.LD.pruned.4DTV.noContig.treebest.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 186_filtered.LD.pruned.4DTV.noContig.treebest.out | awk '{printf $0}' > 186_filtered.LD.pruned.4DTV.noContig.treebest.nwk
treebest nj -b 1000 186_filtered.LD.pruned.noContig.min4.fasta > 186_filtered.LD.pruned.noContig.treebest.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 186_filtered.LD.pruned.noContig.treebest.out | awk '{printf $0}' > 186_filtered.LD.pruned.noContig.treebest.nwk


iqtree2 -s 186_filtered.LD.pruned.4DTV.noContig.min4.fasta  -m MFP --alrt 1000 -B 1000 -T 20 -keep-ident & 
iqtree2 -s 186_filtered.LD.pruned.noContig.min4.fasta       -m MFP --alrt 1000 -B 1000 -T 20 -keep-ident & 

