cd  /root/workspace/186sample/PopLDdecay

PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop East.pop -OutStat East.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop West.pop -OutStat West.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop Central.pop -OutStat Central.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop Central-East.pop -OutStat Central-East.PopLDdecay &

Plot_OnePop.pl -inFile East.PopLDdecay.stat.gz --output East.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West.PopLDdecay.stat.gz --output West.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Central.PopLDdecay.stat.gz --output Central.PopLDdecay.stat -bin1 10 -bin2 100

vi K2_multi.list
West.PopLDdecay.stat.gz West
Central-East.PopLDdecay.stat.gz Central-East

vi K3_multi.list
West.PopLDdecay.stat.gz West
Central.PopLDdecay.stat.gz Central
East.PopLDdecay.stat.gz East


Plot_MultiPop.pl -inList K2_multi.list --output K2_West&Central-East.PopLDdecay -bin1 10 -bin2 100
Plot_MultiPop.pl -inList K3_multi.list --output K3_West&Central&East.PopLDdecay -bin1 10 -bin2 100

