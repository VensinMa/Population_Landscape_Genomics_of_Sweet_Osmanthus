
mkdir /public1/guop/mawx/workspace/cultivar_snpcalling/fastq.gz/raw_data
find /public1/guop/mawx/workspace/cultivar_snpcalling/fastq.gz/ -type f -name "*.fq.gz" -exec cp {} /public1/guop/mawx/workspace/cultivar_snpcalling/fastq.gz/raw_data \;
mv /public1/guop/mawx/workspace/cultivar_snpcalling/fastq.gz/raw_data  /public1/guop/mawx/workspace/cultivar_snpcalling/1.raw_data 
cd /public1/guop/mawx/workspace/cultivar_snpcalling/1.raw_data
ls | awk -F'_' '{print $1}' | sort | uniq | less -SN 




