if(!requireNamespace("lfmm", quietly = TRUE)) {  
  remotes::install_github("bcm-uga/lfmm")
}
if (!requireNamespace("qvalue", quietly = TRUE)) {
  BiocManager::install("qvalue")
}
if (!requireNamespace("LEA", quietly = TRUE)) {
  BiocManager::install("data.table")
}
if (!requireNamespace("LEA", quietly = TRUE)) {
  BiocManager::install("LEA")
}

library(lfmm)
library(qvalue)  
library(data.table)
library(LEA)


# 设置工作目录
getwd()
setwd("C:/Rstudio/RStudio/Workspace/GEA_2024/lfmm")
getwd()

#ped2lfmm("194samples_snp.nounanchor.renamed.filtered.vcftools.plink.ped", 
 #        "194samples_filtered.nounanchor.lfmm", force = TRUE)
######################  读入解释变量与响应变量  ###############################
## X Y 数据中样品顺序（行）要一致
## 读取环境数据 （解释变量）
X = read.csv("/public1/guop/mawx/R/workspace/GEA_2024/lfmm/Climate_current_194samples.csv", 
             header = T, row.names = 1)
X = X[,3:25]
str(X)
## 读取基因型数据（响应变量）sed '1d; s/NA/9/g' 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.recodeA.raw | awk '{ $1=$2=$3=$4=$5=$6=""; print substr($0, index($0,$7)) }' > 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.lfmm
Y = fread("/public1/guop/mawx/R/workspace/GEA_2024/lfmm/194samples_snp.nounanchor.renamed.filtered.vcftools.plink.lfmm")
K = 4  # admixture确定的最佳分群数

mod.lfmm <- lfmm::lfmm_ridge(Y = Y,
                             X = X,
                             K = K)

pv <- lfmm::lfmm_test(Y = Y, X = X, 
                      lfmm = mod.lfmm,
                      calibrate = "gif")

# 保存与加载 pv
saveRDS(pv, file = "pv_results.rds")
pv <- readRDS("pv_results.rds")
##################### 原始p/校正p 值 的提取与保存 ##############################
# 读取校正后的pv$calibrated.pvalue 作为p值
raw_pvalues <- pv$pvalue
calibrated_pvalues <- pv$calibrated.pvalue

pv$gif

# 读取 SNP 标识符文件   awk '{print $1":"$4}' 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.map > 194samples_snp.nounanchor.renamed.filtered.vcftools.plink.ID
snp_ids <- fread("/public1/guop/mawx/R/workspace/GEA_2024/lfmm/194samples_snp.nounanchor.renamed.filtered.vcftools.plink.ID", 
                 header = FALSE)

# 检查 p 值矩阵的维度是否与 SNPid 文件一致
if (nrow(raw_pvalues) != nrow(snp_ids)) {
  stop("行数不匹配：SNPid 文件与 p 值矩阵的行数不同。")
}

# 将 SNPid 作为行名，而不是第一列
rownames(raw_pvalues) <- snp_ids$V1
rownames(calibrated_pvalues) <- snp_ids$V1

# 分别合并 SNPid 和 原始p/校正p 值
result_raw_pvalues <- data.frame(raw_pvalues = raw_pvalues)
result_calibrated_pvalues <- data.frame(calibrated_pvalues = calibrated_pvalues)

# 保存合并后的数据为 CSV 文件，设置 row.names = TRUE 以包含 SNPid 作为行名
fwrite(result_raw_pvalues, file = "snp_raw_pvalues_merged.csv", sep = ",", row.names = TRUE)
fwrite(result_calibrated_pvalues, file = "snp_calibrated_pvalues_merged.csv", sep = ",", row.names = TRUE)


############# 校正后的 p 值进一步校正为 q 值（假发现率，FDR）###################

# 对每个环境因子的校正p值进行q值（FDR）校正
num_env_factors <- ncol(calibrated_pvalues)
qvalues_list <- list()  # 用于存储每个环境因子的 q 值

for (i in 1:num_env_factors) {
  # 对每个环境因子的p值进行q值校正
  pvals <- calibrated_pvalues[, i]
  
  # 计算 q 值
  qobj <- qvalue(p = pvals)
  qvals <- qobj$qvalues
  
  # 将 q 值存储到列表中
  qvalues_list[[i]] <- data.frame(SNPid = snp_ids$V1, qvalues = qvals)
  
  # 保存每个环境因子的 q 值为单独的文件
  fwrite(qvalues_list[[i]], file = paste0("qvalues_env_factor_", i, ".csv"), sep = ",", row.names = FALSE)
}

# 如果想要将所有环境因子的 q 值放到一个数据框中，可以进行如下合并：
all_qvalues_df <- do.call(cbind, lapply(qvalues_list, function(x) x$qvalues))
rownames(all_qvalues_df) <- snp_ids$V1  # 设置行名为 SNPid

# 保存所有环境因子的 q 值为一个综合文件
fwrite(all_qvalues_df, file = "all_qvalues_combined.csv", sep = ",", row.names = TRUE)


###################### FDR < 0.05 和 FDR < 0.01 显著性位点的筛选 ###########################
library(qvalue)
library(data.table)

# 读取 SNP 标识符文件
snp_ids <- fread("C:/Users/mawenxin/Downloads/Linux_down/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink.snp.id", header = FALSE)

# 获取环境因子的列名
env_factor_names <- colnames(X)  # 这是解释变量 X 的列名

# 创建列表保存每个环境因子的阈值、对应的 FDR（q 值），以及显著性位点的数量
threshold_list <- list()
significant_counts <- data.frame(Environment_Factor = env_factor_names, Significant_SNP_Count = integer(num_env_factors))

# 对每个环境因子的校正 p 值进行 q 值（FDR）校正，并找到 q 值等于 0.05 和 0.01 的 p 值阈值
num_env_factors <- ncol(calibrated_pvalues)

for (i in 1:num_env_factors) {
  # 获取每个环境因子的校正 p 值
  pvals <- calibrated_pvalues[, i]
  
  # 计算 q 值
  qobj <- qvalue(p = pvals)
  qvals <- qobj$qvalues
  
  # 找到 q 值最接近 0.05 和 0.01 的校正 p 值
  diff_05 <- abs(qvals - 0.05)
  closest_index_05 <- which.min(diff_05)  # 找到与 0.05 最接近的 q 值的索引
  
  diff_01 <- abs(qvals - 0.01)
  closest_index_01 <- which.min(diff_01)  # 找到与 0.01 最接近的 q 值的索引
  
  # 获取相应的校正 p 值和 q 值作为阈值
  threshold_value_05 <- pvals[closest_index_05]
  corresponding_qvalue_05 <- qvals[closest_index_05]
  
  threshold_value_01 <- pvals[closest_index_01]
  corresponding_qvalue_01 <- qvals[closest_index_01]
  
  # 保存阈值和对应的 FDR（q 值）到列表中，使用环境因子的名称而不是索引
  threshold_list[[i]] <- data.frame(Environment_Factor = env_factor_names[i], 
                                    Threshold_05 = threshold_value_05, 
                                    FDR_qvalue_05 = corresponding_qvalue_05,
                                    Threshold_01 = threshold_value_01,
                                    FDR_qvalue_01 = corresponding_qvalue_01)
  
  # 保存所有位点的 SNPid、p 值和 q 值
  all_results_df <- data.frame(SNPid = snp_ids$V1, 
                               pvalues = pvals, 
                               qvalues = qvals)
  fwrite(all_results_df, file = paste0("all_results_", env_factor_names[i], ".csv"), sep = ",", row.names = FALSE)
  
  # 筛选出 FDR < 0.05 的显著性位点
  significant_indices <- which(qvals < 0.05)
  
  if (length(significant_indices) > 0) {
    significant_pvals <- pvals[significant_indices]
    significant_qvals <- qvals[significant_indices]
    significant_snp_ids <- snp_ids$V1[significant_indices]
    
    # 保存显著性位点的 SNPid、p 值和 q 值
    significant_results_df <- data.frame(SNPid = significant_snp_ids, 
                                         pvalues = significant_pvals, 
                                         qvalues = significant_qvals)
    fwrite(significant_results_df, file = paste0("significant_results_", env_factor_names[i], ".csv"), sep = ",", row.names = FALSE)
    
    # 单独保存显著性位点的 SNPid
    fwrite(as.data.table(significant_snp_ids), file = paste0("significant_snp_ids_", env_factor_names[i], ".csv"), sep = ",", row.names = FALSE)
    
    # 统计 FDR < 0.05 的显著性位点数量
    significant_counts[i, "Significant_SNP_Count"] <- length(significant_snp_ids)
  } else {
    significant_counts[i, "Significant_SNP_Count"] <- 0
  }
}

# 将所有环境因子的阈值和对应的 FDR（q 值）保存为 CSV 文件
threshold_df <- do.call(rbind, threshold_list)
fwrite(threshold_df, file = "environment_factors_thresholds_at_qvalue_0.05_and_0.01.csv", sep = ",", row.names = FALSE)

# 保存显著性位点数量
fwrite(significant_counts, file = "significant_snp_counts.csv", sep = ",", row.names = FALSE)

###################### 统计所有环境因子关联到的位点数量（并集） ###########################

# 创建一个集合来存储所有环境因子下显著性位点的并集
all_significant_snp_ids <- unique(unlist(lapply(1:num_env_factors, function(i) {
  # 获取每个环境因子的校正 p 值
  pvals <- calibrated_pvalues[, i]
  
  # 计算 q 值
  qobj <- qvalue(p = pvals)
  qvals <- qobj$qvalues
  
  # 筛选出 FDR < 0.05 的显著性位点的 SNPid
  significant_indices <- which(qvals < 0.05)
  
  if (length(significant_indices) > 0) {
    return(snp_ids$V1[significant_indices])
  } else {
    return(NULL)
  }
})))

# 将显著性位点的 SNPid 保存为数据框，并包含显著性位点的数量
all_significant_snp_df <- data.frame(SNPid = all_significant_snp_ids)
significant_snp_count <- nrow(all_significant_snp_df)

# 保存显著性位点及其数量为 CSV 文件
fwrite(all_significant_snp_df, file = "all_significant_snp_ids_union.csv", sep = ",", row.names = FALSE)

# 打印显著性位点的数量
cat("所有环境因子下显著性位点的数量（并集）为：", significant_snp_count, "\n")

## 所有环境因子下显著性位点的数量（并集）为： 28603   31381

################################ 绘制曼哈顿图 ##################################

library(CMplot)
library(data.table)

# 读取数据
data <- fread("all_results_BIO19.csv")  # 使用实际的文件名

# 分割 SNPid 列，提取 CHR 和 POSITION
data[, `:=`(CHR = as.numeric(sub("Chr(\\d+):.*", "\\1", SNPid)),  # 提取染色体编号
            POSITION = as.numeric(sub(".*:(\\d+)", "\\1", SNPid)))]  # 提取位置信息

# 创建一个包含所需列的数据框：SNPID, CHR, POSITION, PVALUES
cmplot_data <- data[, .(SNPID = SNPid, CHR, POSITION, PVALUES = pvalues)]

# SNPs 位置和基因名称，外显子和内含子 SNP
exonic_snps <- c("Chr23:13769389", "Chr23:4976511", "Chr11:12632521", "Chr09:23138593",
                 "Chr18:8651829", "Chr11:17984450", "Chr21:24185162", "Chr04:18272634",
                 "Chr03:24194792", "Chr23:9808413", "Chr01:12502472", "Chr09:5697056",
                 "Chr05:24200854", "Chr20:11218738", "Chr15:25977055", "Chr10:22808719",
                 "Chr17:16801696", "Chr15:18779813")
exonic_genes <- c("WRKY9", "JMT", "RBL", "LEA", "RLK6", "IF3", "VAL1", "ENTH", "PEPPER",
                  "TCF21", "TFIIB", "SHT", "PPR", "LPLAT1", "ERF109", "SPL1", "CKI", "SKOR")

# 内含子 SNP
intronic_snps <- c("Chr18:18790486", "Chr17:17815908", "Chr10:6572779", "Chr02:26658889",
                   "Chr15:25270809", "Chr02:33437633", "Chr16:29708694", "Chr06:2071789")
intronic_genes <- c("ARP4", "MSH1", "ATXR2", "PO", "CK2", "SOX2", "BCAT", "UMPS")

# 合并外显子和内含子 SNP 及基因
all_snps <- c(exonic_snps, intronic_snps)
all_genes <- c(exonic_genes, intronic_genes)

# 设置高亮符号，外显子和内含子使用不同符号
highlight_pch <- c(rep(19, length(exonic_snps)), rep(17, length(intronic_snps)))  # 外显子圆形，内含子三角形

# 生成曼哈顿图并在图中展示两条 FDR 阈值线
pdf("PL_plot_BIO19.pdf", width = 15, height = 5)

CMplot(cmplot_data,
       plot.type = "m",
       LOG10 = TRUE,  # 将 P 值进行 log10 转换
       col = c("#423D77", "#3F678B", "#468C8D", "#5FB47F", "#9FD55C", "#F9E956"),  # 点颜色
       
       # 设置 Y 轴范围
       ylim = c(0, 12.5),  # Y 轴显示范围最高到 12
       
       # 高亮所有 SNP（外显子和内含子）
       highlight = all_snps,  # 高亮 SNP 合并后的列表
       highlight.col = c(rep("red", length(exonic_snps)), rep("blue", length(intronic_snps))),  # 外显子 SNP 使用红色，内含子 SNP 使用蓝色
       highlight.cex = 0.7,  # 高亮点大小相同
       highlight.pch = highlight_pch,  # 外显子和内含子使用不同的符号
       
       # 为所有 SNP 添加基因名称标注
       highlight.text = all_genes,
       highlight.text.col = "black",  # 基因名称颜色为黑色
       highlight.text.cex = 1,  # 基因名称文字大小
       
       # 添加 FDR 阈值线
       threshold = c(0.000255424016988622, 0.0000213857275814149),  # 两个 FDR 阈值
       threshold.lty = c(2, 2),  # 两条线的虚线样式
       threshold.col = c("gray20", "gray50"),  # 设置第一条阈值线为黑色，第二条为灰色
       threshold.lwd = c(2, 2),  # 两条线的宽度
       
       # 其他设置
       file.output = FALSE,  # 不输出额外文件
       amplify = FALSE,  # 不放大显著 SNP
       verbose = F,  # 关闭冗长输出
       cex = 0.5  # 点大小
)

dev.off()

#########################  曼哈顿图y轴显示范围调整  ###########################

library(CMplot)
library(data.table)

# 读取数据
data <- fread("all_results_BIO19.csv")  # 使用实际的文件名

# 分割 SNPid 列，提取 CHR 和 POSITION
data[, `:=`(CHR = as.numeric(sub("Chr(\\d+):.*", "\\1", SNPid)),  # 提取染色体编号
            POSITION = as.numeric(sub(".*:(\\d+)", "\\1", SNPid)))]  # 提取位置信息

# 计算 -log10(P-values)
data[, log10P := -log10(pvalues)]

# 过滤 -log10(P-values) > 12.5 的数据点
cmplot_data <- data[log10P <= 12.5, .(SNPID = SNPid, CHR, POSITION, PVALUES = pvalues)]

# SNPs 位置和基因名称，外显子和内含子 SNP
exonic_snps <- c("Chr23:13769389", "Chr23:4976511", "Chr11:12632521", "Chr09:23138593",
                 "Chr18:8651829", "Chr11:17984450", "Chr21:24185162", "Chr04:18272634",
                 "Chr03:24194792", "Chr23:9808413", "Chr01:12502472", "Chr09:5697056",
                 "Chr05:24200854", "Chr20:11218738", "Chr15:25977055", "Chr10:22808719",
                 "Chr17:16801696", "Chr15:18779813")
exonic_genes <- c("WRKY9", "JMT", "RBL", "LEA", "RLK6", "IF3", "VAL1", "ENTH", "PEPPER",
                  "TCF21", "TFIIB", "SHT", "PPR", "LPLAT1", "ERF109", "SPL1", "CKI", "SKOR")

# 内含子 SNP
intronic_snps <- c("Chr18:18790486", "Chr17:17815908", "Chr10:6572779", "Chr02:26658889",
                   "Chr15:25270809", "Chr02:33437633", "Chr16:29708694", "Chr06:2071789")
intronic_genes <- c("ARP4", "MSH1", "ATXR2", "PO", "CK2", "SOX2", "BCAT", "UMPS")

# 合并外显子和内含子 SNP 及基因
all_snps <- c(exonic_snps, intronic_snps)
all_genes <- c(exonic_genes, intronic_genes)

# 设置高亮符号，外显子和内含子使用不同符号
highlight_pch <- c(rep(19, length(exonic_snps)), rep(17, length(intronic_snps)))  # 外显子圆形，内含子三角形

# 生成曼哈顿图并在图中展示两条 FDR 阈值线
pdf("PL_plot_BIO19.pdf", width = 15, height = 5)

CMplot(cmplot_data,
       plot.type = "m",
       LOG10 = TRUE,  # 将 P 值进行 log10 转换
       col = c("#423D77", "#3F678B", "#468C8D", "#5FB47F", "#9FD55C", "#F9E956"),  # 点颜色
       
       # 设置 Y 轴范围
       ylim = c(0, 12.5),  # Y 轴显示范围最高到 12.5
       
       # 高亮所有 SNP（外显子和内含子）
       highlight = all_snps,  # 高亮 SNP 合并后的列表
       highlight.col = c(rep("red", length(exonic_snps)), rep("blue", length(intronic_snps))),  # 外显子 SNP 使用红色，内含子 SNP 使用蓝色
       highlight.cex = 0.7,  # 高亮点大小相同
       highlight.pch = highlight_pch,  # 外显子和内含子使用不同的符号
       
       # 为所有 SNP 添加基因名称标注
       highlight.text = all_genes,
       highlight.text.col = "black",  # 基因名称颜色为黑色
       highlight.text.cex = 1,  # 基因名称文字大小
       
       # 添加 FDR 阈值线
       threshold = c(0.000255424016988622, 0.0000213857275814149),  # 两个 FDR 阈值
       threshold.lty = c(2, 2),  # 两条线的虚线样式
       threshold.col = c("gray20", "gray80"),  # 设置第一条阈值线为黑色，第二条为灰色
       threshold.lwd = c(2, 2),  # 两条线的宽度
       
       # 其他设置
       file.output = FALSE,  # 不输出额外文件
       amplify = FALSE,  # 不放大显著 SNP
       verbose = F,  # 关闭冗长输出
       cex = 0.5  # 点大小
)

dev.off()

##########################################################################
library(data.table)

# 读取 significant_results_BIO19.csv 文件
significant_results <- fread("C:/Rstudio/RStudio/Workspace/GEA/lfmm/significant_results_BIO5.csv")

# 读取 1861SNP_NR注释副本.csv 文件，假设第三列是包含 SNP 的列
annotation_data <- fread("E:/OneDrive/文档/1861SNP_NR注释副本.csv")

# 重命名 1861SNP_NR注释副本.csv 的第三列为 "SNPid"，以便后续交集操作
setnames(annotation_data, old = names(annotation_data)[3], new = "SNPid")

# 提取与 significant_results 中 SNPid 列有交集的行
matched_significant_results <- significant_results[SNPid %in% annotation_data$SNPid]
matched_annotation_data <- annotation_data[SNPid %in% significant_results$SNPid]

# 合并两个数据集以便保存
combined_data <- merge(matched_significant_results, matched_annotation_data, by = "SNPid")

# 检查数据结构
head(combined_data)

# 使用 fwrite 保存并指定编码为 UTF-8
fwrite(combined_data, file = "E:/OneDrive/文档/matched_significant_results_BIO5.csv", 
       sep = ",", 
       bom = TRUE,  # 添加 BOM 标记，确保 UTF-8 被正确识别
       quote = FALSE)
