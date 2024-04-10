cd  /root/workspace/186sample/PopLDdecay

PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop East.pop -OutStat East.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop West.pop -OutStat West.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop Central.pop -OutStat Central.PopLDdecay &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -SubPop Central-East.pop -OutStat CentralEast.PopLDdecay &

Plot_OnePop.pl -inFile East.PopLDdecay.stat.gz --output East.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West.PopLDdecay.stat.gz --output West.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Central.PopLDdecay.stat.gz --output Central.PopLDdecay.stat -bin1 10 -bin2 100

vi K2_multi.list
West.PopLDdecay.stat.gz West
CentralEast.PopLDdecay.stat.gz CentralEast

vi K3_multi.list
West.PopLDdecay.stat.gz West
Central.PopLDdecay.stat.gz Central
East.PopLDdecay.stat.gz East


Plot_MultiPop.pl -inList K2_multi.list --output K2_West2CentralEast.PopLDdecay -bin1 10 -bin2 100
Plot_MultiPop.pl -inList K3_multi.list --output K3_West2Central2East.PopLDdecay -bin1 10 -bin2 100


####################################################################################################################################

cd  /root/workspace/186sample/PopLDdecay

PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -MaxDist 10 -SubPop East.pop -OutStat East.PopLDdecay_10K &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -MaxDist 10 -SubPop West.pop -OutStat West.PopLDdecay_10K &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -MaxDist 10 -SubPop Central.pop -OutStat Central.PopLDdecay_10K &
PopLDdecay  -InVCF /mnt/e/mwx/workspace/pop/186_filtered_vcftools.noContig.recode.vcf -MaxDist 10 -SubPop Central-East.pop -OutStat CentralEast.PopLDdecay_10K &

Plot_OnePop.pl -inFile East.PopLDdecay_10K.stat.gz --output East.PopLDdecay_10K.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West.PopLDdecay_10K.stat.gz --output West.PopLDdecay_10K.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Central.PopLDdecay_10K.stat.gz --output Central.PopLDdecay_10K.stat -bin1 10 -bin2 100

vi K2_multi_10K.list
West.PopLDdecay_10K.stat.gz West
CentralEast.PopLDdecay_10K.stat.gz CentralEast

vi K3_multi_10K.list
West.PopLDdecay_10K.stat.gz West
Central.PopLDdecay_10K.stat.gz Central
East.PopLDdecay_10K.stat.gz East


Plot_MultiPop.pl -inList K2_multi_10K.list --output K2_West2CentralEast.PopLDdecay_10K -bin1 10 -bin2 100
Plot_MultiPop.pl -inList K3_multi_10K.list --output K3_West2Central2East.PopLDdecay_10K -bin1 10 -bin2 100
