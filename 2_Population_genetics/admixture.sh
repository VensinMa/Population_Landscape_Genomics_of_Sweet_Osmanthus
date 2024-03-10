# 需要bed格式文件
plink --vcf 224__filtered.LD.pruned.noContig.recode.vcf --make-bed   --out 224_filtered.LD.pruned.noContig --keep-allele-order --allow-extra-chr  

# 计算K=2到20
seq 2 20 | parallel -j 20 "admixture --cv  224_filtered.LD.pruned.noContig.bed {} 1>admix.{}.log 2>&1"

# 确定具有最小CV值的K值为最佳分群数
cat *.log | grep "CV"

mkdir result
cp  ./*.Q result/

# 绘制admixture遗传结构一键脚本 可在linux windows环境运行 需注意软件、文件路径
Rscript  draw_admixture.R  result  224.filtered.LD.pruned.noContig.nosex   224.filtered.LD.pruned.noContig.admixture
Rscript.exe  draw_admixture.R result 224.filtered.LD.pruned.noContig.nosex  224.filtered.LD.pruned.noContig.admixture


######### 推荐本地运行 ###### Rstudio中 ########## 便于调色校正 #################
# 安装需要的包并加载
if (!requireNamespace("argparser", quietly = TRUE)) install.packages("argparser")
if (!requireNamespace("pophelper", quietly = TRUE)) install.packages("pophelper")
if (!requireNamespace("cols4all", quietly = TRUE)) install.packages("cols4all")
library(pophelper)
library(argparser)
library(cols4all)

# 设置工作路径 输入文件路径
setwd("C:/RStudio/RStudio/Workspace/Resequence/admixture/")
dir <- "result" # Q矩阵文件所在的目录
sample_file <- "224.filtered.LD.pruned.noContig.nosex" # 样本顺序文件路径
output_prefix <- "224.filtered.LD.pruned.noContig" # 输出文件前缀

# 定义颜色调色板
c4a_gui()
set.seed(12)
mycol = c4a("classic_green_orange12", 12)
# 随机打乱颜色向量
mycol1 <- sample(mycol)
mycol1 = c("#55ACD9","#D4994D","#DC5166","#F2B8E8","#86DCD3","#6C559D",
           "#B7CE76","#FFE7A9","#F0AF75","#B1D5BA","#E68A8B","#D6CDE0")

# 读取Q矩阵文件
Qfiles <- list.files(dir, pattern="Q", full.names=TRUE)
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



# 单独绘制某一k值
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

