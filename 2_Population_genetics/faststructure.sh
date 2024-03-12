##  fastStructure使用的输入文件基本与admixture 一致，bed格式

cd /public1/guop/mawx/workspace/pop_genetic/faststructure
mkdir fast_result
seq 2 30 | parallel -j 30 "structure.py -K {} --input=/public1/guop/mawx/workspace/pop_genetic/224.filtered.LD.pruned.noContig --output=fast_result/224.filtered.LD.pruned.noContig_{} --cv=5 --prior=logistic --seed=100 > log.{}.txt 2>&1" &

###################################### 本地绘图 ######################################
# 安装需要的包并加载
if (!requireNamespace("argparser", quietly = TRUE)) install.packages("argparser")
if (!requireNamespace("pophelper", quietly = TRUE)) install.packages("pophelper")
if (!requireNamespace("cols4all", quietly = TRUE)) install.packages("cols4all")
library(pophelper)
library(argparser)
library(cols4all)

# 设置工作路径 输入文件路径
setwd("C:/RStudio/RStudio/Workspace/Resequence/faststructure")
dir <- "faststructureRes/" # Q矩阵文件所在的目录
sample_file <- "wild_fragrans_filtered_LD.nosex" # 样本顺序文件路径
output_prefix <- "wild_fragrans_filtered_LD.nosex" # 输出文件前缀

# 定义颜色调色板
c4a_gui()
set.seed(12)
mycol = c4a("palette36", 36)
# 随机打乱颜色向量
mycol1 <- sample(mycol)

# 读取Q矩阵文件
Qfiles <- list.files(dir, pattern="meanQ", full.names=TRUE)
qlist <- readQ(Qfiles, indlabfromfile=FALSE)

# 读取样本标签
sample_labels <- read.table(sample_file, header=FALSE)
for (i in 1:length(qlist)) {
  rownames(qlist[[i]]) <- sample_labels$V1
}

# 绘制所有k值结构
plotQ(qlist,
      sortind="all", # 对所有个体进行排序
      imgtype="pdf", # 图像类型，可以是"pdf", "png"等
      ordergrp=FALSE, # 是否按群组排序
      imgoutput="join", # 图像输出方式，"join"为合并所有Q矩阵的图形
      width=40, # 图像宽度
      height=8, # 图像高度
      titlecol = "black",
      clustercol=mycol1,
      exportpath=getwd(), # 图像导出路径，这里设置为当前工作目录
      outputfilename=output_prefix, # 输出文件名前缀
      useindlab=TRUE, # 是否使用个体标签
      sharedindlab=FALSE, # 标签是否共享
      showindlab=TRUE) # 是否显示个体标签



# 单独绘制某一k值 直接导出
plotQ(qlist["224.filtered.LD.pruned.noContig.6.Q"],
      sortind="all", # 对所有个体进行排序
      imgtype="pdf", # 图像类型，可以是"pdf", "png"等
      ordergrp=FALSE, # 是否按群组排序
      width=40, # 图像宽度
      height=8, # 图像高度
      titlecol = "black",
      clustercol=mycol1,
      exportpath=getwd(), # 图像导出路径，这里设置为当前工作目录
      outputfilename=output_prefix, # 输出文件名前缀
      useindlab=TRUE, # 是否使用个体标签
      sharedindlab=FALSE, # 标签是否共享
      showindlab=TRUE) # 是否显示个体标签

# 单独绘制某一k值 在r中
plotQ(qlist["wild_fragrans_filtered_LD_output.6"],
      sortind="all", # 对所有个体进行排序
      clustercol=mycol1,
      ordergrp=F, # 是否按群组排序
      width=40, # 图像宽度
      height=8, # 图像高度
      titlecol = "black",
      exportplot=FALSE,
      returnplot=TRUE,
      useindlab=TRUE, # 是否使用个体标签
      sharedindlab=FALSE, # 标签是否共享
      showindlab=TRUE) # 是否显示个体标签


# 单独绘制某一k值 直接导出
plotQ(qlist,
      sortind="all",
      imgtype="pdf",
      ordergrp=FALSE,
      imgoutput="join",
      width=40,
      height=8,
      exportpath=getwd(),
      outputfilename="test0",
      useindlab=TRUE,
      sharedindlab=FALSE,
      showindlab=TRUE)

