library(psych) 


# 设置工作目录
getwd()
setwd("C:/Rstudio/RStudio/Workspace/GEA_2024/lfmm")
getwd()

## 加载群体经纬度环境数据
env_data = read.csv("extracted_data/Climate_current_194samples.csv",
                    header = T, row.names = 1)

## 相关性绘图
pairs.panels(env_data[, 3:25], smooth = TRUE, scale = TRUE, density = TRUE,
             method = "pearson", pch = 20, lm = FALSE, cor = TRUE,
             jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE,
             rug = TRUE, breaks = "Sturges",  wt = NULL, ellipses = TRUE, 
             cex.cor = 5, cex = 1, smoother = FALSE, stars = TRUE, alpha = .05, 
             hist.border = "black")
xx = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
  "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

pairs.panels(env_data[, xx], 
             smooth = TRUE, scale = TRUE, density = TRUE,
             method = "pearson", pch = 20, lm = FALSE, cor = TRUE,
             jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE,
             rug = TRUE, breaks = "Sturges",  wt = NULL, ellipses = TRUE, 
             cex.cor = 5, cex = 0.8, smoother = FALSE, stars = TRUE, alpha = .05, 
             hist.border = "black")

########################################
# 加载必要的包
library(psych)

# 读取数据文件
env_data <- read.csv("extracted_data/Climate_current_194samples.csv", header = TRUE, row.names = 1)

# 计算皮尔逊相关系数矩阵
cor_matrix <- cor(env_data[, 3:25], method = "pearson", use = "pairwise.complete.obs")

# 获取相关性绝对值大于0.8的变量对
high_cor <- which(abs(cor_matrix) > 0.8, arr.ind = TRUE)

# 过滤掉对角线元素以及重复的组合
high_cor <- high_cor[high_cor[,1] != high_cor[,2], ]
high_cor <- high_cor[high_cor[,1] < high_cor[,2], ]

# 打印相关性绝对值大于0.8的变量对及其相关系数
if (nrow(high_cor) > 0) {
  for (i in 1:nrow(high_cor)) {
    var1 <- rownames(cor_matrix)[high_cor[i, 1]]
    var2 <- colnames(cor_matrix)[high_cor[i, 2]]
    cor_value <- cor_matrix[high_cor[i, 1], high_cor[i, 2]]
    cat(sprintf("Variables: %s and %s, Correlation: %.2f\n", var1, var2, cor_value))
  }
} else {
  cat("No variable pairs have a correlation above 0.8.\n")
}
