mkdir -p /home/vensin/workspace/PCA
cd /home/vensin/workspace/PCA

## vcf转bed格式
plink --vcf /home/vensin/workspace/population_structure/admixture/194samples_filtered.LD.pruned.Superscaffold2Chr.recode.vcf \
    --make-bed   --out /home/vensin/workspace/PCA/194samples_filtered.LD.pruned.Superscaffold2Chr  --keep-allele-order  --allow-extra-chr  

## 生成grm文件
gcta64 --bfile 194samples_filtered.LD.pruned.Superscaffold2Chr  --autosome  --make-grm  --out GA

## 进行PCA分析
gcta64 --grm GA --pca 130  --out PCA_out

## 一键绘图脚本 draw_PCA.R 可在2_Population_genetics目录获得
Rscript ./draw_PCA.R  PCA_out.eigenvec 1 2  sample.pop  PCA_out_figure
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

# 输入参数设置
input <- "PCA_out.eigenvec" # eigenvec文件路径
x <- 1 # X轴对应的主成分编号
y <- 2 # Y轴对应的主成分编号
pop_file <- "sample.pop" # 群体信息文件路径
outpre <- "PCA_out_figure" # 输出文件前缀

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
pdf(file=paste(outpre, "pc", x, y, "pdf", sep = "."), height = 10, width = 12)
print(p1)
dev.off()

png(file=paste(outpre, "pc", x, y, "png", sep = "."), height = 3000, width = 3600)
print(p1)
dev.off()
