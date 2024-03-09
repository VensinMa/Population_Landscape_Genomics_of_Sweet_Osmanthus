args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量是否正确
if(length(args) != 5) {
  cat("Usage: Rscript draw_PCA.R <PCA_out.eigenvec> <PC1_index> <PC2_index> <sample.pop> <output_prefix>\n")
  cat("Example: Rscript draw_PCA.R PCA_out.eigenvec 1 2 sample.pop PCA_out_figure\n")
  quit(status = 1) # 退出状态1表示有错误发生
}

# 参数赋值
input <- args[1]
x <- as.numeric(args[2])
y <- as.numeric(args[3])
pop_file <- args[4]
outpre <- args[5]

# 检查并安装缺失的包
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}
if (!requireNamespace("ggrepel", quietly = TRUE)) {
  install.packages("ggrepel")
}
# if (!requireNamespace("ggsci", quietly = TRUE)) {
#   install.packages("ggsci")
# }

library(ggplot2)
library(ggrepel)
# library(ggsci)

# 创建输出文件夹
output_dir <- "PCA_draw_out"
if(!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# 读取数据
vec <- read.table(input, header = FALSE, row.names = 1, sep = " ")
pop <- read.table(pop_file, header = FALSE, row.names = 1, sep = "\t")

# 准备数据
population <- pop[rownames(vec), 1]
individual <- pop[rownames(vec), 2]

vec$population <- population
vec$individual <- individual

# 绘图
p1 <- ggplot(data = vec, aes(x = vec[,x+1], y = vec[,y+1], colour = population)) +
  geom_point() +
  #geom_text_repel(aes(label = individual), size = 5, box.padding = 0.35, point.padding = 0.5, max.overlaps = Inf) +
  #geom_text(aes(label = individual), vjust = 1.5, hjust = 0.5, size = 3) +  # 添加标签
  theme_bw() +
  xlab(paste("PC", x, sep = "")) +
  ylab(paste("PC", y, sep = ""))

# 构造完整的PDF和PNG文件路径
pdf_file_path <- paste(output_dir, paste(outpre, "pc", x, y, "pdf", sep = "."), sep = "/")
png_file_path <- paste(output_dir, paste(outpre, "pc", x, y, "png", sep = "."), sep = "/")

# 保存图形为PDF
pdf(pdf_file_path, height = 10, width = 12)
print(p1)
dev.off()

# 保存图形为PNG
png(png_file_path, height = 3000, width = 3600, res = 500)
print(p1)
dev.off()



