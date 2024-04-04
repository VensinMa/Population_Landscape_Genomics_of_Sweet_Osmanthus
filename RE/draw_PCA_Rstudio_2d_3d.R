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
library(dplyr)
library(plotly)
library(scatterplot3d)
library(RColorBrewer)
library(cols4all)

setwd("C:/RStudio/RStudio/Workspace/Resequence/186sample/PCA")
# 输入参数设置
input <- "186_filtered.LD.pruned.noContig_plink_pca.eigenvec" # eigenvec文件路径
x <- 1 # X轴对应的主成分编号
y <- 2 # Y轴对应的主成分编号
z <- 3

pop_file <- "186sampleid_K2.txt" # 群体信息文件路径
outpre <- "186_LD_plink_PCA_out_figure_K2" # 输出文件前缀

# 读取数据
vec <- read.table(input, header = T, row.names = 1, sep = " ")
pop <- read.table(pop_file, header = FALSE, row.names = 1, sep = "\t")

# 准备数据
population <- pop[rownames(vec), 1]
individual <- pop[rownames(vec), 2]

vec$population <- population
vec$individual <- individual

# 读取特征值文件
eigenvalues <- read.table("186_filtered.LD.pruned.noContig_plink_pca.eigenval")

# 计算特征值总和
total_variance <- sum(eigenvalues$V1)

# 计算每个主成分的占比
proportions <- eigenvalues$V1 / total_variance

# 打印占比
print(proportions)

pc1_variance <- "56.0%"  # PC1解释的方差百分比
pc2_variance <- "24.1%"   # PC2解释的方差百分比
pc3_variance <- "19.9%"   # PC3解释的方差百分比

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
p1 <- ggplot(data = vec, aes(x = vec[,x+1], y = vec[,y+1], colour = population)) +
  geom_point() +
  scale_color_manual(values = c("#FFD700","#7FB2D5")) +  # 使用ggsci的Nature颜色方案
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


# 使用plotly创建三维散点图
p <- plot_ly(data = vec, x = ~vec[,x+1], y = ~vec[,y+1], z = ~vec[,z+1], color = ~population, 
             text = ~individual, marker = list(size = 5)) %>%
  add_markers() %>%
  layout(scene = list(xaxis = list(title = paste("PC", x)),
                      yaxis = list(title = paste("PC", y)),
                      zaxis = list(title = paste("PC", z))),
         title = "PCA 3D Plot")
p
marker = list(size = 2)
library(plotly)

# 假设 vec 是您的数据框，其中包含 x+1、y+1 和 z+1 列作为主成分分析的结果，
# population 是颜色分组的变量，individual 是每个点的标签。

p <- plot_ly(data = vec, 
             x = ~vec[,x+1], 
             y = ~vec[,y+1], 
             z = ~vec[,z+1], 
             color = ~population, 
             text = ~individual,
             type = "scatter3d", 
             mode = "markers",
             marker = list(size = 3)) %>%
  layout(scene = list(
    xaxis = list(title = paste("PC", x, pc1_variance), showgrid = TRUE, gridcolor = 'black', zeroline = TRUE),
    yaxis = list(title = paste("PC", y, pc2_variance), showgrid = TRUE, gridcolor = 'black', zeroline = TRUE),
    zaxis = list(title = paste("PC", z, pc3_variance), showgrid = TRUE, gridcolor = 'black', zeroline = TRUE)
  ),
  title = "PCA 3D Plot with Grid")

p
str(p)



##################################       #####################################

# 这里我们使用Set1调色板，您可以根据需要选择其他的调色板
# colors <- brewer.pal(n = length(unique(vec$population)), name = "Set1")

c4a_gui()
vec$population <- as.factor(vec$population)
set.seed(11)
colors = c4a("dynamic", length(unique(vec$population)))
# 使用 sample() 函数随机打乱 colors 中的颜色
shuffled_colors <- sample(colors)
colors =  shuffled_colors
# 将因子类型的population列映射到颜色向量
color_vector <- colors[as.numeric(vec$population)]

# 绘制散点图
s3d <- scatterplot3d(vec[, x+1], vec[, y+1], vec[, z+1], 
                     color = color_vector,  # 使用映射后的颜色向量
                     pch = 16,
                     angle=30,
                     cex.symbols = 1,
                     type = "p",
                     main = "PCA 3D Plot with Grid",
                     xlab = paste("PC", x, "(", pc1_variance, "%)", sep=""),
                     ylab = paste("PC", y, "(", pc2_variance, "%)", sep=""),
                     zlab = paste("PC", z, "(", pc3_variance, "%)", sep=""),
                     grid = TRUE,
                     box = T)

s3d
