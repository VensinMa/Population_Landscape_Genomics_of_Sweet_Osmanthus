#!/bin/bash

# 基本参数设置
INPUT_FOLDER="/path/to/input"        # 输入文件夹
OUTPUT_FOLDER="/path/to/output"     # 输出文件夹
GENOME_FOLDER="/path/to/genome"     # 基因组文件夹
REFERENCE="${GENOME_FOLDER}/bwa/hg19.fa"
REF_FLAT="${GENOME_FOLDER}/refFlat_new.txt"
ACCESS_BED="${GENOME_FOLDER}/access-5k-mappable.hg19.bed"

# 线程数
THREADS=30

# 日志文件
LOG_FILE="${OUTPUT_FOLDER}/process.log"
mkdir -p $OUTPUT_FOLDER
echo "Process started at $(date)" > $LOG_FILE

# 函数：记录日志和输出到终端
log_and_echo() {
  echo "$1" | tee -a $LOG_FILE
}

# 遍历所有的 *.qc.r1.fastq.gz 文件，自动匹配对应的 R2 文件
for R1 in ${INPUT_FOLDER}/*.qc.r1.fastq.gz; do
  # 获取样本名（去掉后缀部分）
  BASENAME=$(basename $R1 .qc.r1.fastq.gz)
  R2="${INPUT_FOLDER}/${BASENAME}.qc.r2.fastq.gz"

  # 检查 R2 文件是否存在
  if [ ! -f $R2 ]; then
    log_and_echo "Error: Missing R2 file for $BASENAME. Skipping..."
    continue
  fi

  # 输出文件前缀
  PREFIX="${OUTPUT_FOLDER}/${BASENAME}"

  log_and_echo "Processing sample: $BASENAME"

  # 质控
  log_and_echo "Step 1: Quality control using fastp for $BASENAME..."
  fastp -i $R1 -I $R2 \
        -o ${PREFIX}.qc.r1.fastq.gz -O ${PREFIX}.qc.r2.fastq.gz \
        --trim_front1 7 --trim_tail1 0 --trim_front2 7 --trim_tail2 0 \
        -h ${PREFIX}.qc.html -j ${PREFIX}.qc.json &>> $LOG_FILE
  log_and_echo "Step 1 completed for $BASENAME."

  # 比对
  log_and_echo "Step 2: Alignment using BWA for $BASENAME..."
  bwa mem -t $THREADS $REFERENCE ${PREFIX}.qc.r1.fastq.gz ${PREFIX}.qc.r2.fastq.gz \
    | samtools view -bS - > ${PREFIX}.bam &>> $LOG_FILE
  log_and_echo "Step 2 completed for $BASENAME."

  # 排序
  log_and_echo "Step 3: Sorting BAM file for $BASENAME..."
  samtools sort ${PREFIX}.bam -o ${PREFIX}.sorted.bam -@ $THREADS &>> $LOG_FILE
  log_and_echo "Step 3 completed for $BASENAME."

  # 索引
  log_and_echo "Step 4: Indexing BAM file for $BASENAME..."
  samtools index ${PREFIX}.sorted.bam -@ $THREADS &>> $LOG_FILE
  log_and_echo "Step 4 completed for $BASENAME."

  # 移除重复
  log_and_echo "Step 5: Removing duplicates using Sambamba for $BASENAME..."
  sambamba markdup -t $THREADS ${PREFIX}.sorted.bam ${PREFIX}.sorted.rm.bam -r &>> $LOG_FILE
  samtools index ${PREFIX}.sorted.rm.bam -@ $THREADS &>> $LOG_FILE
  log_and_echo "Step 5 completed for $BASENAME."

  # CNV 分析
  log_and_echo "Step 6: Performing CNV analysis for $BASENAME..."
  cnvkit.py batch ${PREFIX}.sorted.rm.bam -n \
      -m wgs --segment-method cbs --drop-low-coverage --target-avg-size 500000 \
      -f $REFERENCE --annotate $REF_FLAT --access $ACCESS_BED \
      --output-reference ${PREFIX}_reference.cnn \
      --output-dir $OUTPUT_FOLDER &>> $LOG_FILE
  log_and_echo "Step 6 completed for $BASENAME."

  # 绘制图表
  log_and_echo "Step 7: Generating CNV scatter and heatmap for $BASENAME..."
  cnvkit.py scatter ${PREFIX}.sorted.rm.cnr -s ${PREFIX}.sorted.rm.cns -o ${PREFIX}_scatter.pdf --y-max 2 --y-min -2 &>> $LOG_FILE
  cnvkit.py heatmap ${PREFIX}.sorted.rm.cns -o ${PREFIX}_heatmap.pdf &>> $LOG_FILE
  log_and_echo "Step 7 completed for $BASENAME."

  # 重复标记与统计
  log_and_echo "Step 8: Marking duplicates and generating statistics for $BASENAME..."
  gatk MarkDuplicates -I ${PREFIX}.sorted.bam -O ${PREFIX}.markup.bam -M ${PREFIX}.metrics.txt &>> $LOG_FILE
  samtools flagstat ${PREFIX}.markup.bam -O ${PREFIX}.flagstat.txt &>> $LOG_FILE
  log_and_echo "Step 8 completed for $BASENAME."

  log_and_echo "Finished processing: $BASENAME"
done

log_and_echo "All samples processed successfully. Process finished at $(date)."
