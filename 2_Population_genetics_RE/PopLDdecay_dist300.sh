cd  /root/workspace/186sample/PopLDdecay

### Step1: 计算 LD decay
PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop East.pop -OutStat East.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop Central.pop -OutStat Central.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop West-GZ.pop -OutStat West-GZ.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop West-YN.pop -OutStat West-YN.PopLDdecay

PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop Central-East.pop -OutStat Central-East.PopLDdecay
PopLDdecay  -InVCF /home/vensin/workspace/PopLDdecay/194_vcftools.recode.vcf -MaxDist 300 -SubPop West.pop -OutStat West.PopLDdecay

### Step2: 绘图
## 单个谱系绘图
Plot_OnePop.pl -inFile East.PopLDdecay.stat.gz --output East.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile Central.PopLDdecay.stat.gz --output Central.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West-GZ.PopLDdecay.stat.gz --output West-GZ.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West-YN.PopLDdecay.stat.gz --output West-YN.PopLDdecay.stat -bin1 10 -bin2 100

Plot_OnePop.pl -inFile Central-East.PopLDdecay.stat.gz --output Central-East.PopLDdecay.stat -bin1 10 -bin2 100
Plot_OnePop.pl -inFile West.PopLDdecay.stat.gz --output West.PopLDdecay.stat -bin1 10 -bin2 100

## 多谱系绘图
cat <<EOF > K4_multi.list
East.PopLDdecay.stat.gz East
Central.PopLDdecay.stat.gz Central
West-GZ.PopLDdecay.stat.gz West-GZ
West-YN.PopLDdecay.stat.gz West-YN
EOF

cat <<EOF > K2_multi.list
Central-East.stat Central-East
West.PopLDdecay.stat West
EOF

Plot_MultiPop.pl -inList K4_multi.list --output K4.PopLDdecay -bin1 10 -bin2 100
Plot_MultiPop.pl -inList K2_multi.list --output K2.PopLDdecay -bin1 10 -bin2 100

gzip -d -k Central-East.PopLDdecay.stat.bin.gz Central.PopLDdecay.stat.bin.gz  East.PopLDdecay.stat.bin.gz  West-GZ.PopLDdecay.stat.bin.gz  West-YN.PopLDdecay.stat.bin.gz  West.PopLDdecay.stat.bin.gz
