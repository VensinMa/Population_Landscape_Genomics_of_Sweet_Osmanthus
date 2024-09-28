

cd 
python /root/software/vcf2phylip-2.8/vcf2phylip.py --input 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.recode.vcf --fasta --output-prefix 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned
python /root/software/vcf2phylip-2.8/vcf2phylip.py --input 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.recode.vcf --fasta --output-prefix 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing




treebest nj -b 1000 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.min4.fasta > 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.treebest.out
sed 's/:\([0-9.]\+\)\[&&NHX:B=\([0-9]\+\)\]/\2:\1/' 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.treebest.out | awk '{printf $0}' > 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.nomissing.treebest.nwk


treebest nj -b 1000 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.min4.fasta > 200samples_filtered_renamed.snp.unanchor.final.vcftools.LD.pruned.treebest.out




