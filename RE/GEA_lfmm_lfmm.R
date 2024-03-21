if(!requireNamespace("lfmm", quietly = TRUE)) {  
  remotes::install_github("bcm-uga/lfmm")
}
if (!requireNamespace("qvalue", quietly = TRUE)) {
  install.packages("qvalue")
}

library(lfmm)
library(qvalue)  
library(data.table)

# 设置工作目录
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
## X Y 数据中样品顺序（行）要一致
## 读取环境数据 （解释变量）
X = read.csv("186sample_id_lat_lon_ele_with_env_vars.csv", 
             header = T, row.names = 1)
X = X[,3:25]
str(X)
## 读取基因型数据（响应变量）
Y = fread("186_filtered.LD.pruned.noContig.recode.lfmm")
K = 2  # admixture确定的最佳分群数

mod.lfmm <- lfmm::lfmm_ridge(Y = Y,
                             X = X,
                             K = K)

pv <- lfmm::lfmm_test(Y = Y, X = X, 
                          lfmm = mod.lfmm,
                          calibrate = "gif")

# 读取校正后的pv$calibrated.pvalue 作为p值
pvalues <- pv$calibrated.pvalue
pv$gif
# 将 pvalues 矩阵保存为 CSV 文件
# write.csv(pvalues, file = "calibrated_pvalues.csv", row.names = FALSE)
fwrite(pvalues, file = "calibrated_pvalues.csv",sep = "," )

## 将校正后的 p 值进一步校正为 q 值（假发现率，FDR）
# pv$calibrated.pvalue 是包含校正后 p 值的矩阵
qvalues_list <- lapply(1:ncol(pv$calibrated.pvalue), function(i) {
  q_vals <- qvalue(pv$calibrated.pvalue[, i])$qvalues
  return(q_vals)
})
colnames(pv$calibrated.pvalue)
# 为整个列表设置名称
names(qvalues_list) <- colnames(pv$calibrated.pvalue)

# 转换列表为矩阵
qvalues_matrix <- do.call(cbind, qvalues_list)

# 设置矩阵列名
colnames(qvalues_matrix) <- names(qvalues_list)
colnames(qvalues_matrix)

# 初始化一个列表来保存每列 q 值小于 0.05 的位点索引
significant_indices_list <- list()

# 遍历 qvalues_matrix 的每一列
for(i in 1:ncol(qvalues_matrix)) {
  # 找到每列中 q 值小于 0.05 的位点的索引
  significant_indices <- which(qvalues_matrix[, i] < 0.05)
  # 将找到的索引保存到列表中
  significant_indices_list[[colnames(qvalues_matrix)[i]]] <- significant_indices
}

# significant_indices_list 现在包含了每列中显著位点的索引

# 计算每列中显著位点的数量
significant_counts <- sapply(significant_indices_list, length)
significant_counts
# significant_counts 是一个向量，包含每列显著位点的数量
'''
    > significant_counts
         BIO1      BIO2      BIO3      BIO4      BIO5      BIO6      BIO7      BIO8 
         1357      1075      9271      1162      2767        93       462       896 
         BIO9     BIO10     BIO11     BIO12     BIO13     BIO14     BIO15     BIO16 
         2312       565        66      2688      1234       525       287       298 
        BIO17     BIO18     BIO19 Elevation      SRAD       SOC     PHH2O 
         1331      7029      9142      2186       864       932      4475 
'''



#####  写出位点名称
# 指定文件夹名称
significant_indices_folder_name <- "significant_indices"
# 检查文件夹是否存在，如果不存在则创建
if (!dir.exists(significant_indices_folder_name)) {
  dir.create(significant_indices_folder_name)
}

# 遍历 significant_indices_list 中的每一列
for (col_name in names(significant_indices_list)) {
  # 生成文件名，包含文件夹路径和列名
  filename <- file.path(significant_indices_folder_name, paste0(col_name, "_significant_indices.txt"))
  
  # 获取当前列的显著位点索引
  significant_indices <- significant_indices_list[[col_name]]
  
  # 将索引保存到文件中
  writeLines(as.character(significant_indices), filename)
}

# 求显著位点索引的并集
all_significant_indices <- unique(unlist(significant_indices_list))
length(all_significant_indices)
'''
    > significant_counts
         BIO1      BIO2      BIO3      BIO4      BIO5      BIO6      BIO7      BIO8 
         1357      1075      9271      1162      2767        93       462       896 
         BIO9     BIO10     BIO11     BIO12     BIO13     BIO14     BIO15     BIO16 
         2312       565        66      2688      1234       525       287       298 
        BIO17     BIO18     BIO19 Elevation      SRAD       SOC     PHH2O 
         1331      7029      9142      2186       864       932      4475 
'''

full_path <- file.path(significant_indices_folder_name, "all_significant_indices.txt")

# 将并集保存到文件中
writeLines(as.character(all_significant_indices), full_path)

#########################  FDR 0.05 的 p 值阈值  ##############################

# 初始化一个向量来保存每列的 p 值阈值
pvalue_thresholds <- numeric(ncol(qvalues_matrix))

# 遍历每列
for(i in 1:ncol(qvalues_matrix)) {
  # 获取当前列
  qvalues <- qvalues_matrix[, i]
  # 找出 q 值 <= 0.05 的最大 q 值
  max_qvalue_under_threshold <- max(qvalues[qvalues <= 0.05])
  
  # 找到与这个 q 值相对应的最小 p 值
  # 注意：这里假设 p 值与 q 值是单调关联的，这在大多数情况下是成立的
  corresponding_pvalue <- min(pv$calibrated.pvalue[qvalues == max_qvalue_under_threshold, i])
  
  # 保存这个 p 值到阈值向量中
  pvalue_thresholds[i] <- corresponding_pvalue
}

pvalue_thresholds
# pvalue_thresholds 向量现在包含了每列在 FDR 0.05 下的 p 值阈值
'''
    > pvalue_thresholds
     [1] 3.771611e-05 2.990105e-05 2.589362e-04 3.234397e-05 7.727829e-05 2.540473e-06
     [7] 1.288894e-05 2.494263e-05 6.551284e-05 1.574773e-05 1.729086e-06 7.503425e-05
    [13] 3.428627e-05 1.461574e-05 7.694276e-06 8.298022e-06 3.715242e-05 1.962802e-04
    [19] 2.553796e-04 6.080596e-05 2.411788e-05 2.594347e-05 1.249710e-04
    > 
'''


# 初始化向量来保存每列的 q 值阈值
qvalue_thresholds <- numeric(ncol(qvalues_matrix))

# 遍历每列，这个循环现在也计算和保存 q 值阈值
for(i in 1:ncol(qvalues_matrix)) {
  qvalues <- qvalues_matrix[, i]
  # 找出 q 值 <= 0.05 的最大 q 值，作为 q 值阈值
  max_qvalue_under_threshold <- max(qvalues[qvalues <= 0.05])
  qvalue_thresholds[i] <- max_qvalue_under_threshold  # 保存 q 值阈值
  
  # 计算对应的 p 值阈值
  corresponding_pvalue <- min(pv$calibrated.pvalue[qvalues == max_qvalue_under_threshold, i])
  pvalue_thresholds[i] <- corresponding_pvalue
}

# 创建一个数据框来保存 p 值和 q 值阈值
thresholds_df <- data.frame(qvalue_thresholds = qvalue_thresholds,
                            pvalue_thresholds = pvalue_thresholds)

# 将阈值数据框写入 CSV 文件
write.csv(thresholds_df, file = "thresholds.csv", row.names = FALSE)




length(which(qvalues_matrix[,2] < 0.10)) 
length(which(qvalues_matrix[,2] < 0.05)) 

FDR.05 <- which(qvalues_matrix[,2] < 0.05)


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

plot(-log10(pv$calibrated.pvalue[,1]), 
     pch = 19, 
     cex = .3,
     xlab = "Probe", 
     ylab = "-Log P",
     col = "grey",
     main = "Column 1")

# 在图中添加对应 FDR 0.05 的阈值横线
threshold_line <- -log10(pvalue_thresholds[1])
abline(h = threshold_line, col = "red", lwd = 2)





plot(-log10(pv$calibrated.pvalue[,18]), 
     pch = 19, 
     cex = .3,
     xlab = "Probe", ylab = "-Log P",
     col = "grey")

points(causal.set, 
       -log10(pv$calibrated.pvalue)[causal.set], 
       col = "blue")

##  3、qq图
qqplot(rexp(length(pvalues), rate = log(10)),
       -log10(pvalues), xlab = "Expected quantile",
       pch = 19, cex = .4)
abline(0,1)






