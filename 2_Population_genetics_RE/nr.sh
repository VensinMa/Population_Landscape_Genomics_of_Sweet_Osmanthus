## NR注释记录 20240325
### 1、安装 aspera 与下载 NR 蛋白数据库文件

# 下载安装 aspera 软件
conda create -n aspera
conda activate aspera
conda install -c hcc aspera-cli -y

# 添加并立即生效环境变量
echo 'export PATH="/home/vensin/anaconda3/envs/aspera/bin/:$PATH"' >> ~/.bashrc
source ~/.bashrc

# 检查是否安装成功
ascp -h

# 切换到工作目录 workspace 
cd /home/vensin/workspace/nr.annotations || mkdir -p /home/vensin/workspace/nr.annotations && cd /home/vensin/workspace/nr.annotations
    
# 下载 nr 蛋白库
ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz ./ &
ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/FASTA/nr.gz.md5 ./ &
# wget https://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz.md5
md5sum -c nr.gz.md5
# gunzip -c nr.gz > nr.fasta
unpigz -c -p 16 nr.gz > nr.fasta



    
# 下载 Nr 数据库中蛋白编号 mprot.accession 和物种编号 taxid 的对应关系信息 
# ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/pub/taxonomy/accession2taxid/prot.accession2taxid.gz ./ &
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz

# 下载 Nr 数据库中物种编号 taxdmp 的层次信息 
# ascp -i /home/vensin/anaconda3/envs/aspera/etc/asperaweb_id_dsa.openssh -l 1000M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
# wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip

cd /home/vensin/data/genome
perl gff_filter_longest.pl LYG.hic.gff longest.cds.ids LYG.longest.cds.gff
gffread  LYG.hic.gff -g LYG.hic.fasta -w LYG.gene.fasta
gffread  LYG.hic.gff  -g  LYG.hic.fasta  -y  LYG.pep.fasta 
gffread  LYG.longest.cds.gff  -g  LYG.hic.fasta  -x  LYG.longest.cds.fasta 
gffread  LYG.longest.cds.gff  -g  LYG.hic.fasta  -y  LYG.longest.pep.fasta 

## (base) vensin@ubuntu24-04:~/data/genome$ awk '{print $3}' LYG.hic.gff | grep -c "gene"
## 41252
## (base) vensin@ubuntu24-04:~/data/genome$ grep -c ">" LYG.gene.fasta
## 56392
## (base) vensin@ubuntu24-04:~/data/genome$ grep -c ">" LYG.pep.fasta
## 56392
## (base) vensin@ubuntu24-04:~/data/genome$ grep -c ">" LYG.longest.pep.fasta
## 41252
## (base) vensin@ubuntu24-04:~/data/genome$ grep -c ">" LYG.longest.cds.fasta
## 41252

### 2、使用 Blast/Diamond 进行 NR 注释

cd /home/vensin/workspace/nr.annotations || mkdir -p /home/vensin/workspace/nr.annotations && cd /home/vensin/workspace/nr.annotations
# 使用diamond  软件的子命令 makedb 将 fasta 格式的蛋白序列创建成后缀为 dmnd 的数据库文件：
diamond makedb --in nr.fasta --db nr.db
# Writing trailer...  [108.497s]
# Closing the input file...  [0.786s]
# Closing the database file...  [0.361s]

# Database sequences  707028945
#   Database letters  272881947790
#      Database hash  3cfcab0838efcbd44c935e5519724c33
#         Total time  8928s

## sudo fallocate -l 50G /swapfile
## sudo chmod 600 /swapfile
## sudo mkswap /swapfile
## sudo swapon /swapfile

# 将物种全基因组核酸/蛋白序列 blastx / blastp 到构建好的数据库：
diamond blastp --db /home/vensin/workspace/nr.annotations/nr.db.dmnd --query /home/vensin/data/genome/LYG.longest.pep.fasta \
    --out /home/vensin/workspace/nr.annotations/LYG.longest.pep.Nr.annotations \
    --outfmt 6 qseqid sseqid pident evalue bitscore qlen slen length mismatch gapopen qstart qend sstart send stitle \
    --sensitive --max-target-seqs 5 --evalue 1e-5 --index-chunks 1 
    # --index-chunks 1 
    
diamond blastp --db /public1/guop/liujj/liujj/workspace/Annovar/nr/nr.db.dmnd \
    --query /public1/guop/liujj/liujj/workspace/Annovar/nr/genome/LYG.longest.pep.fasta \
    --out   \
    --outfmt 6 qseqid sseqid pident evalue bitscore qlen slen length mismatch gapopen qstart qend sstart send stitle \
    --sensitive --max-target-seqs 5 --evalue 1e-5 --index-chunks 1 --threads 64

diamond blastp --db /public1/guop/liujj/liujj/workspace/Annovar/nr/nr.db.dmnd \
    --query /public1/guop/liujj/liujj/workspace/Annovar/nr/genome/LYG.longest.pep.fasta \
    --out /public1/guop/liujj/liujj/workspace/Annovar/nr/genome/LYG.longest.pep.Nr.annotations \
    --outfmt 6 qseqid sseqid pident evalue bitscore qlen slen length mismatch gapopen qstart qend sstart send stitle \
    --sensitive --max-target-seqs 5 --evalue 1e-5 --index-chunks 1 --threads 64
## Total time = 69305.8s
## Reported 198346 pairwise alignments, 198346 HSPs.
## 40044 queries aligned.

# diamond 默认设置下输出表格格式的结果。结果分12列，其结果信息和 BLAST 默认设置-outfmt 6输出的格式完全一致。
#  1. qseqid     query序列ID 
#  2. sseqid     subject序列ID
#  3. pident     Identity百分比
#  4. length     匹配长度
#  5. mismatch   错配长度	
#  6. gapopen    打开Gap的次数
#  7. qstart     query序列起始位点
#  8. qend       query序列结束位点
#  9. sstart     subject序列起始位点
# 10. send       subject序列结束位点
# 11. evalue     E-vaule值
# 12. bitscore   bitscore得分
    
'''    
# 使用 blast 软件的命令 makeblastdb 将 fasta 格式的蛋白序列创建成后缀为.psp .pin .phr 的数据库文件：
makeblastdb -in nr.test.fasta -dbtype prot -out nr.test
    
# 将物种全基因组核酸/蛋白序列 blastx / blastp 到构建好的数据库：    
blastx -query Of.tps.test.fasta -db nr.test -evalue 1e-5 -max_target_seqs 1 -outfmt '6 qseqid sseqid pident evalue bitscore qlen slen length mismatch gapopen qstart qend sstart send stitle' -out Of.tps.test

# blastp用蛋白质序列搜索蛋白质序列库;blastn用核酸序列搜索核酸库;blastx核酸序列对蛋白质库的比对，核酸序列在比对之前自动按照六个读码框翻译成蛋白质序列。
# -db 指定blast搜索用的数据库，后接我们自己建的库，或者下载解压好的NR，NT等数据库。
# -evalue 1e-5 设置e值。
# -max_target_seqs 1 设置最多的目标序列匹配数，相当于可视化blast比对结果中，红色的序列个数。
# -outfmt format 6  5代表xml格式，这个格式的信息最全，6代表table格式，该格式的结果以table的形式给出，清晰易懂，7代表有注释行的table格式，比6多了一些#开头的注释行。还有一种可以自己定义输出内容，例如：在输出6格式的前提下，输出一列为物种名字： -outfmt "6 std ssciname"。
 '''   
    
    
