if(!requireNamespace("lfmm", quietly = TRUE)) {  
  remotes::install_github("bcm-uga/lfmm")
}
if (!requireNamespace("qvalue", quietly = TRUE)) {
  install.packages("qvalue")
}

library(lfmm)
library(qvalue)  
library(data.table)

## 参考教程链接 https://bookdown.org/hhwagner1/LandGenCourse_book/WE_11.html
# 目标目录
getwd()
dir_lfmm <- "lfmm"

# 检查目录是否存在
if (!file.exists(dir_lfmm)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_lfmm)
}
setwd(dir_lfmm)
getwd()
######################  读入解释变量与响应变量  ###############################
## 读取环境数据 （解释变量）
X = read.csv("186sample_id_lat_lon_ele_with_env_vars.csv", 
             header = T, row.names = 1)
X = X[,3:25]

## 读取基因型数据（响应变量）
Y = fread("186_filtered.LD.pruned.noContig.recode.lfmm")
K = 5  # admixture确定的最佳分群数

mod.lfmm <- lfmm::lfmm_ridge(Y = Y,
                             X = X,
                             K = K)

pv <- lfmm::lfmm_test(Y = Y, X = X, 
                          lfmm = mod.lfmm,
                          calibrate = "gif")

# 读取校正后的pv$calibrated.pvalue 作为p值
pvalues <- pv$calibrated.pvalue
pv$gif

## 将校正后的 p 值进一步校正为 q 值（假发现率，FDR）
# 使用 lapply 转换 p 值为 q 值，并保留列名
qvalues_list <- lapply(1:ncol(pv$calibrated.pvalue), function(i) {
  col_name <- colnames(pv$calibrated.pvalue)[i] # 获取当前列名
  q_vals <- qvalue(pv$calibrated.pvalue[, i])$qvalues # 计算q值
  setNames(list(q_vals), col_name) # 将q值列表赋予列名
})
# 将列表转换为矩阵
qvalues_matrix <- do.call(cbind, qvalues_list)
# 由于列表的每个元素名已经是原始列名，我们可以直接将这些名字设置为矩阵的列名
colnames(qvalues_matrix) <- names(qvalues_list)


length(which(qv < 0.10)) ## h.w many SNPs have an FDR < 10%?
length(which(qv < 0.05)) ## h.w many SNPs have an FDR < 5%?
(FDR.1 <- colnames(Y)[which(qv < 0.1)]) ## i.entify which SNPs these are
(FDR.05 <- colnames(Y)[which(qv < 0.05)])

write.csv(FDR.1, "FDR_0.1_SNPs.csv", quote = FALSE, row.names = FALSE)
write.csv(FDR.05, "FDR_0.05_SNPs.csv", quote = FALSE, row.names = FALSE)


################################  绘图  #####################################

## 1、p值分布直方图

# 设置图形布局为 2 行 1 列
par(mfrow = c(2, 1))
## 未校正 p值分布直方图 
hist(pv$pvalue[,18], main="Unadjusted p-values")   
## 校正后 p值分布直方图 
hist(pv$calibrated.pvalue[,18], main="GIF-adjusted p-values")
# 重置图形参数到默认值
par(mfrow = c(1, 1))

# 遍历每一列
for (i in 1:ncol(pv$pvalue)) {
  # 使用列名作为文件名
  col_name <- colnames(pv$pvalue)[i]
  pdf_filename <- paste(col_name, "_pvalue_distributions.pdf", sep = "")
  
  # 开始将图形输出到 PDF 文件
  pdf(pdf_filename)
  
  # 设置图形布局为 2 行 1 列
  par(mfrow = c(2, 1))
  
  # 绘制未校正 p 值分布直方图，不显示 X 轴标题
  hist(pv$pvalue[, i], main = paste("Unadjusted p-values -", col_name), xlab = "")   
  
  # 绘制校正后 p 值分布直方图，不显示 X 轴标题
  hist(pv$calibrated.pvalue[, i], main = paste("GIF-adjusted p-values -", col_name), xlab = "")
  
  # 重置图形参数到默认值
  par(mfrow = c(1, 1))
  
  # 关闭 PDF 文件
  dev.off()
}


##  2、曼哈顿图
plot(-log10(pv$calibrated.pvalue[,18]), 
     pch = 19, 
     cex = .3,
     xlab = "Probe", ylab = "-Log P",
     col = "grey")
causal.set <- seq(11, 1496, by = 80)
points(causal.set, 
       -log10(pv$calibrated.pvalue)[causal.set], 
       col = "blue")

##  3、qq图
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)






