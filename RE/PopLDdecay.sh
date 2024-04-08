

PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop east.pop -OutStat east.PopLDdecay
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop west.pop -OutStat west.PopLDdecay

Plot_OnePop.pl -inFile east.PopLDdecay.stat.gz --output east.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile west.PopLDdecay.stat.gz --output west.PopLDdecay.stat -bin1 10 -bin2 100

vi multi.list
east.PopLDdecay.stat.gz east
west.PopLDdecay.stat.gz west

Plot_MultiPop.pl -inList multi.list --output east_west.PopLDdecay -bin1 10 -bin2 100
