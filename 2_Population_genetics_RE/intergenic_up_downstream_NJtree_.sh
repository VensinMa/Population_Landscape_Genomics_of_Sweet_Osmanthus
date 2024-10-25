mkdir -p /home/vensin/workspace/intergenic_up_downstream_NJtree
cd /home/vensin/workspace/intergenic_up_downstream_NJtree
cp /home/vensin/workspace/Annovar/194samples_filtered.intergenic_up_downstream.LD.prune.recode.vcf /home/vensin/workspace/intergenic_up_downstream_NJtree
VCF2Dis -i  194samples_filtered.intergenic_up_downstream.LD.prune.recode.vcf  -o 194samples_filtered.intergenic_up_downstream.LD.prune_dis.mat    -Rand 0.25 &

