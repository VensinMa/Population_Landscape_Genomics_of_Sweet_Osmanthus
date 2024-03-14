fasttree -nt  229_filtered.LD.pruned.noContig.min4.fasta > 229_filtered.LD.pruned.noContig.fasttree.out


iqtree -s 229_filtered.LD.pruned.noContig.min4.fasta  -m MFP  -B 1000  -keep-ident  -T 20 &
