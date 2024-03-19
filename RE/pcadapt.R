
library(pcadapt)
library(LEA)
library(cols4all)

getwd()
# 目标目录
dir_name <- "pcadapt"
# 检查目录是否存在
if (!file.exists(dir_name)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_name)
}
# 改变工作目录到目标目录
setwd(dir_name)
getwd()

# vcf2lfmm("Z:/root/workspace/186sample/186_filtered.LD.pruned.noContig.recode.vcf", "186_filtered.LD.pruned.noContig.lfmm", force = TRUE)

#########################  1. 读取基因型数据 ###################################
path_to_file <- "186_filtered.LD.pruned.noContig.recode.lfmm"
# path_to_file <- system.file("extdata", "geno3pops.bed", package = "pcadapt")
filename <- read.pcadapt(path_to_file, type = "lfmm")

#####################   2. 选择主成分的数量 K  ################################
# 用足够大的主成分数量（例如 K=20）首先进行主成分分析
x <- pcadapt(input = filename, K = 20) 
# x <- pcadapt(input = filename, K = 3) 

# 2.1. Scree Plot（碎石图）
# "scree plot" 按降序显示每个 PC 解释的方差百分比
plot(x, option = "screeplot")

# 2.2. Score Plot（得分图）
# 另一种选择 K 的方法是基于显示种群结构的“得分图”
# 得分图显示了个体在指定主成分上的投影

# With names
poplist.names <- c(rep("CP", 6), rep("CX", 5), rep("DA", 6), rep("DRS", 6), 
                   rep("DST", 7), rep("EJ", 6), rep("FZA", 6), rep("GHX", 7),
                   rep("HJLL", 5), rep("HYX", 6), rep("JD", 6), rep("JMX", 3), 
                   rep("LCJ", 6), rep("LS", 7), rep("LX", 6), rep("LQ", 2), 
                   rep("QDH", 10), rep("RX", 6), rep("SFZ", 6), rep("SK", 6),
                   rep("SLZ", 6), rep("SL", 5), rep("ST", 8), rep("SXK", 6), 
                   rep("WYL", 4), rep("XC", 6), rep("LQ", 3), rep("XNF", 6),
                   rep("YK", 5), rep("YX", 5), rep("YZY", 6), rep("ZJP", 6), 
                   rep("ZJS", 2))

print(poplist.names)
#c4a_gui()
mycol = c4a("palette36", 32)
mycol

plot(x, option = "scores", pop = poplist.names)
#plot(x, option = "scores", pop = poplist.names, col = mycol)
## 最佳K值计算分组
x <- pcadapt(filename, K = 5)
summary(x)

## 曼哈顿图显示了 p 值的 -log10
plot(x, option = "manhattan")  ## 曼哈顿图
## 使用 Q-Q 图检查 p 值的预期均匀分布
plot(x, option = "qqplot")      ## QQ图
## 测试统计量和 p 值的直方图
hist(x$pvalues, xlab = "p-values", main = NULL, breaks = 50, col = "orange") # p 值的直方图
plot(x, option = "stat.distribution") # 测试统计量的直方图

##################  选择异常检测的截止值  三种方法 ##############################
# q-values 方法
library(qvalue)
qval <- qvalue(x$pvalues)$qvalues
alpha <- 0.05
outliers <- which(qval < alpha)
length(outliers)
hist(qval, xlab = "p-values", main = NULL, breaks = 50, col = "orange") # p 值的直方图

# Benjamini-Hochberg 程序  保守
padj <- p.adjust(x$pvalues, method = "BH",)
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
# Bonferroni 校正  更保守
padj <- p.adjust(x$pvalues, method = "bonferroni")
alpha <- 0.05
outliers <- which(padj < alpha)
length(outliers)
##  > length(outliers)
##  [1] 27412
hist(padj, xlab = "p-values", main = NULL, breaks = 50, col = "orange") # p 值的直方图
outliers
write.csv(outliers, file = "pcadapt_fdr0.05_outliers.csv")




