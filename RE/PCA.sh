## vcf转bed格式
plink --vcf 186_filtered.LD.pruned.noContig.recode.vcf  --make-bed   --out 186_filtered.LD.pruned.noContig  --keep-allele-order  --allow-extra-chr  
plink --vcf 186.snpEffAnno.4dtv.LD.vcf  --make-bed   --out 186_filtered.LD.pruned.4DTV.noContig  --keep-allele-order  --allow-extra-chr 
## 生成grm文件
gcta64 --bfile 186_filtered.LD.pruned.noContig  --autosome  --make-grm  --out GA
gcta64 --bfile 186_filtered.LD.pruned.4DTV.noContig --autosome  --make-grm  --out GA_4D
## 进行PCA分析
gcta64 --grm GA --pca 30  --out 186_filtered.LD.pruned.noContig.PCA_out
gcta64 --grm GA_4D --pca 30  --out 186_filtered.LD.pruned.4DTV.noContig.PCA_out
## 一键绘图脚本 draw_PCA.R 可在2_Population_genetics目录获得
## "Usage: Rscript draw_PCA.R <PCA_out.eigenvec> <PC1_index> <PC2_index> <sample.pop> <output_prefix>"
## "Example: Rscript draw_PCA.R PCA_out.eigenvec 1 2 sample.pop PCA_out_figure"

#########  或选择本地Rstudio绘制 ############### 以下为完整代码
#############################################################

# 安装及加载所需R包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
library(ggsci)
library(ggplot2)
library(ggrepel)


setwd("C:/RStudio/RStudio/Workspace/Resequence/186sample/")
# 输入参数设置
input <- "Z:/root/workspace/186sample/186_filtered.LD.pruned.noContig.PCA_out.eigenvec" # eigenvec文件路径
x <- 1 # X轴对应的主成分编号
y <- 2 # Y轴对应的主成分编号
pop_file <- "186sampleid_pop.txt" # 群体信息文件路径
outpre <- "186_LD_PCA_out_figure" # 输出文件前缀

# 读取数据
vec <- read.table(input, header = FALSE, row.names = 1, sep = " ")
pop <- read.table(pop_file, header = FALSE, row.names = 1, sep = "\t")

# 准备数据
population <- pop[rownames(vec), 1]
individual <- pop[rownames(vec), 2]

vec$population <- population
vec$individual <- individual

p1 <- ggplot(data = vec, aes(x = vec[,x+1], y = vec[,y+1], colour = population)) +
  geom_point() +
  #scale_color_aaas() +  # 使用ggsci的Nature颜色方案
  geom_text_repel(aes(label = individual),  # 使用geom_text_repel来防止文本重叠
                  size = 5,                # 设置文本大小
                  box.padding = 0.35,      # 设置文本周围的空间
                  point.padding = 0.5,     # 设置文本和点之间的空间
                  #segment.color = 'black'
                  max.overlaps = Inf) +  # 设置连接线的颜色
  # geom_text(aes(label = individual), vjust = 1.5, hjust = 0.5, size = 3) +  # 添加标签
  # stat_ellipse() +
  theme_bw() +
  xlab(paste("PC", x, sep = "")) +
  ylab(paste("PC", y, sep = ""))

# 绘制图形
print(p1)
# 保存图形为PDF和PNG
pdf(file=paste(outpre, "pc", x, y, "pdf", sep = "."), height = 20, width = 24)
print(p1)
dev.off()

png(file=paste(outpre, "pc", x, y, "png", sep = "."), height = 3000, width = 3600)
print(p1)
dev.off()
