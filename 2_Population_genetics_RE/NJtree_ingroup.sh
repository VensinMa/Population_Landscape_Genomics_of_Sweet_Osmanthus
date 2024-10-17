cd mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/  && mkdir NJtree
cp 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf ./NJtree
cp 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.log ./NJtree



cd mawx/workspace/wild_snpcalling/4.gatk_gvcf/merged_vcf/194sample/NJtree
python /public1/guop/mawx/software/vcf2phylip-2.8/vcf2phylip.py --input 194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.recode.vcf --fasta --output-prefix 194samples_filtered.LD.pruned
treebest nj -b 1000 194samples_filtered.LD.pruned.min4.fasta > 194samples_filtered.LD.pruned.treebest.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 194samples_filtered.LD.pruned.treebest.out | awk '{printf $0}' > 194samples_filtered.LD.pruned.treebest.nwk

