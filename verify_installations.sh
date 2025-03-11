#!/bin/bash

# 定义需要验证的软件列表
softwares=(
    fastqc bowtie2 multiqc fastp bwa samtools picard bedtools snpeff bcftools gatk
    vcftools blast mafft seqtk diamond hmmer meme iqtree phylip fasttree treebest
    trimal seqkit emboss gffread plink faststructure muscle PopLDdecay xpclr aspera
)

# 遍历软件列表并验证
for software in "${softwares[@]}"; do
    echo "验证 $software..."
    # 激活环境
    conda activate "$software"
    # 检查软件是否可执行
    if command -v "$software" &> /dev/null; then
        echo "$software 安装成功！"
        # 运行帮助命令以进一步验证
        "$software" --help || "$software" -h || echo "$software 无帮助命令"
    else
        echo "$software 未找到，可能安装失败！"
    fi
    # 退出当前环境
    conda deactivate
    echo "-----------------------------"
done
