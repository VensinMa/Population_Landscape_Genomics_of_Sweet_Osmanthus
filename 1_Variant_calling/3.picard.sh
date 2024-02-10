 # 运行Picard MarkDuplicates 注意picard.jar路径
    java -jar /public1/guop/mawx/software/picard/picard.jar MarkDuplicates -I "$output_dir/$species.q30.sorted.bam" -O "$output_dir/$species.markdup.bam" -M "$output_dir/$species.markdup.bam.mat" -MAX_FILE_HANDLES 1000 --REMOVE_DUPLICATES false --TMP_DIR tmp

    # 运行samtools index
    samtools index "$output_dir/$species.markdup.bam"

    # 运行GATK HaplotypeCaller
    gatk --java-options '-Xmx32g -XX:ParallelGCThreads=16 -Djava.io.tmpdir=./tmp' HaplotypeCaller -R "$reference_genome" -I "$output_dir/$species.markdup.bam" -O "$gatk_dir/$species_raw.gvcf" --native-pair-hmm-threads 16 -ERC GVCF --dont-use-soft-clipped-bases

