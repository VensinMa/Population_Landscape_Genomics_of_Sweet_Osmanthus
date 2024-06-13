
cd /root/workspace/stacks

###   (base) root@DESKTOP-GH94RAH:~/workspace/stacks# populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_32pop.txt  --smooth -t 8
###   Error: Malformed arguments: input mode 'vcf' requires an output directory (--out-path).

mkdir populations_output
populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_32pop.txt  --smooth -t 8  -O ./populations_output




cd  /public1/guop/mawx/workspace/186sample
mkdir populations_output
populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_32pop.txt  --smooth -t 40  -O ./populations_output

cd  /public1/guop/mawx/workspace/186sample
mkdir populations_k2_output
populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_K2.txt  --smooth -t 40  -O ./populations_k2_output
mkdir populations_k3_output
populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_K3.txt  --smooth -t 40  -O ./populations_k3_output
mkdir populations_k1_output
populations -V 186_filtered.LD.pruned.noContig.recode.vcf -M 186sample_K1.txt  --smooth -t 40  -O ./populations_k1_output
