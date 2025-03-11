#!/bin/bash

# 初始化 Conda
eval "$(conda shell.bash hook)"

# 更新conda到最新版本，不提示，自动确认
conda update -n base -c defaults conda -y
conda config --add channels defaults

# 定义需要安装的软件列表
softwares=(
    fastqc bowtie2 multiqc fastp bwa samtools picard bedtools snpeff bcftools gatk
    vcftools blast mafft seqtk diamond hmmer meme iqtree phylip fasttree treebest
    trimal seqkit emboss gffread plink faststructure muscle PopLDdecay xpclr aspera
)

# 获取 Conda 基础路径
conda_base=$(conda info --base)

# 遍历软件列表并为每个软件创建一个同名环境
for software in "${softwares[@]}"; do
    # 检查是否已经存在同名环境
    if conda info --envs | grep -q "^${software}\s"; then
        echo "环境 $software 已经存在，跳过..."
    else
        echo "创建环境 $software 并安装相应软件..."
        # 创建一个以软件名称为环境名称的conda环境
        if ! conda create -n "$software" -y; then
            echo "创建环境 $software 失败，跳过..."
            continue
        fi
        # 激活新创建的环境
        conda activate "$software"
        # 安装该软件
        if ! conda install -c bioconda "$software" -y; then
            echo "安装 $software 失败，跳过..."
            conda deactivate
            continue
        fi
        # 退出当前环境
        conda deactivate
    fi

    # 将当前软件环境的路径加入到PATH变量中
    echo "export PATH=\"$conda_base/envs/$software/bin:\$PATH\"" >> ~/.bashrc
done

# 提示用户重新加载bashrc以应用更改
echo "请运行 'source ~/.bashrc' 来应用更改，或重启您的终端。"
