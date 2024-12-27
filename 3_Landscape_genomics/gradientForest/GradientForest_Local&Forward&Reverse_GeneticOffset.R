# 加载所需的R包
# install.packages("conformal", repos="http://R-Forge.R-project.org")
# install.packages("extendedForest", repos="http://R-Forge.R-project.org")
# install.packages("gradientForest", repos="http://R-Forge.R-project.org")
# devtools::install_github("AndiKur4/MaizePal")
library(gradientForest)
library(MaizePal)
library(data.table)
library(gdm)
library(dplyr)
library(tidyverse)
require(raster)
library(fields)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)
library(ggplot2)

setwd("C:/Rstudio/RStudio/Workspace/gradientForest_2024")
#setwd("/public1/guop/mawx/workspace/R/gradientForest_2024")
getwd() # [1] "C:/Rstudio/RStudio/Workspace/gradientForest_2024"
# 设置输入文件、输出结果和图片路径
result_dir <- "final_result"
picture_dir <- "final_picture"
Local_GO_dir <- "Local_Genetic_Offset"
Forward_GO_dir <- "Forward_Genetic_Offset"
Reverse_GO_dir <- "Reverse_Genetic_Offset"
if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(picture_dir)) dir.create(picture_dir)
if (!dir.exists(Local_GO_dir)) dir.create(Local_GO_dir)
if (!dir.exists(Forward_GO_dir)) dir.create(Forward_GO_dir)
if (!dir.exists(Reverse_GO_dir)) dir.create(Reverse_GO_dir)

#############################  GradientForest 模型的构建 #########################################

# 导入群体适应性位点等位基因频率数据
PopsMaf <- fread(file = "input/GF_PopsMaf.csv")
setnames(PopsMaf, old = "V1", new = "pop") # 第一列是群体名称但没有列名 将默认列名V1更改为pop
setDF(PopsMaf) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(PopsMaf) <- PopsMaf[[1]] # 将第一列作为行名
PopsMaf <- PopsMaf[,-1] # 移除已经设置为行名的列
PopsMaf <- PopsMaf[order(rownames(PopsMaf)), ]
PopsMaf <- PopsMaf %>%
  select(where(~ !any(is.na(.))))  # 去除缺失位点
dim(PopsMaf)

# 选择最终保留的用于构建GF模型的预测环境因子
PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

# 读取群体当前环境数据
CurrentEnvs <- read.csv('input/32pop_means_env_vars.csv', 
                        header = T, row.names = 1)
CurrentEnvs <- CurrentEnvs[, PredictEnvs]
PredictEnvsNames = colnames(CurrentEnvs)
print(PredictEnvsNames)

# 计算maxLevel，用于树的最大深度
maxLevel <- log2(0.368 * nrow(PopsMaf) / 2)
Envs_Maf = cbind(CurrentEnvs[, PredictEnvsNames], PopsMaf)

# 构建梯度森林模型，将环境数据和等位基因频率合并在一起
gf.mod <- gradientForest(Envs_Maf,
                         predictor.vars = colnames(CurrentEnvs),# 指定用作预测变量X的列名，即环境数据的列名
                         response.vars = colnames(PopsMaf),# 指定用作响应变量Y的列名，即等位基因频率数据的列名
                         ntree = 1000, maxLevel = maxLevel, trace = T,
                         corr.threshold = 0.50,  nbin = 1001, check.names = FALSE)

# 可以将保存 GF模型结果 gf.mod 保存到到一个文件 便于后续直接加载使用
# save(gf.mod, file = "gf.mod.9674.RData")
# load(file = "gf.mod.9674.RData")

########################################### 绘图 ###################################################

# 绘制预测变量重要性排序图
#生成空的PDF文件 #生成重要值排序 #保存生成的结果
pdf(file = paste0(picture_dir, "/gf.mod.Importance.pdf"), width = 8, height = 8) 
plot(gf.mod, plot.type = "Overall.Importance", 
     col = c(rep("grey",8), MaizePal::maize_pal("HighlandMAGIC", 3)),
     las = 2, cex.names = 0.8) 
dev.off() 

# 保存梯度森林模型的一些输出结果
# 包括响应变量 Y、预测变量 X、重要性值 imp.rsq、模型结果 result等
write.table(gf.mod$Y, file = paste0(result_dir, "/gf.mod.Y.txt"))
write.table(gf.mod$X, file = paste0(result_dir, "/gf.mod.X.txt"))
write.table(gf.mod$imp.rsq, file = paste0(result_dir, "/gf.mod.imp.rsq.txt"))
write.table(gf.mod$result, file = paste0(result_dir, "/gf.mod.result.txt"))
write.table(gf.mod$res.u, file = paste0(result_dir, "/gf.mod.res.u.txt"))
write.table(gf.mod$res, file = paste0(result_dir, "/gf.mod.res.txt"))

# splits density plots 分割密度图plot.gradientForest
PredictEnvs # [1] "BIO2"  "BIO8"  "BIO9"  "BIO10" "BIO12" "BIO15" "BIO17" "BIO18" "SRAD"  "SOC"   "PHH2O"
pdf(file = paste0(picture_dir, "/splits.density.plots_BIO17.pdf"), 
    width = 10, height = 7)
plot(gf.mod, plot.type= "S", imp.vars = "BIO18", leg.posn = "topright", 
     cex.legend = 1.5, cex.axis = 1.2, cex.lab = 1.5, line.ylab = 0.2, 
     par.args = list(mgp = c(2, 0.5, 0), mar = c(3,2,0.1,0.5)))
dev.off()
summary(gf.mod$X$BIO17)
sum(is.na(gf.mod$X$BIO17))

hist(gf.mod$X$BIO17, main="Histogram of BIO17", xlab="BIO17", breaks=20)
# 循环遍历 PredictEnvs 中的每个变量
for (env_var in PredictEnvs) {
  pdf(file = paste0(picture_dir, "/splits.density.plots_", env_var, ".pdf"), 
      width = 10, height = 7)
  plot(gf.mod, plot.type = "S", imp.vars = env_var, leg.posn = "topright", 
       cex.legend = 1.5, cex.axis = 1.2, cex.lab = 1.5, line.ylab = 0.2, 
       par.args = list(mgp = c(2, 0.5, 0), mar = c(3, 2, 0.1, 0.5)))
  dev.off()
}

# 全部环境因子  图例只标注前5个最具响应性的snp
pdf(file = paste0(picture_dir, "/All_env_factor_cumulativeplot_common.scale_F.pdf"), 
    width = 16, height = 18)
plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = "BIO17", 
     show.overall = T, legend = T, leg.posn = "topleft", leg.nspecies = 5,
     common.scale = F, cex.lab = 1.4, cex.legend = 1.2, cex.axis = 1.4, 
     line.ylab = 1.8, par.args = list(mgp = c(2, 0.8, 0), mar = c(4, 2.2, 0.2, 1),
                                      omi = c(0, 0.6, 0.1, 0)))
# mgp 距轴的距离=c(轴标签，轴刻度标签，轴刻度)；内边距 mar=c(下，左，上，右)；外边距 omi
dev.off()

## 单个环境因子
pdf(file = paste0(picture_dir, "/PredictEnvs_cumulativeplot_SRAD.pdf"), 
    width = 5, height = 7)
plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = "SRAD", leg.nspecies = 5,
     show.overall = T, legend = T, leg.posn = "topleft", common.scale = T, 
     cex.lab = 1, cex.legend = 0.7, cex.axis = 1, line.ylab = 2,
     par.args = list(mgp = c(2, 0.5, 0), mar = c(0.5, 0.2, 0.2, 1),
                     omi = c(0.3, 0.7, 0.1, 0.1)))
# mgp 距轴的距离=c(轴标签，轴刻度标签，轴刻度)；内边距 mar=c(下，左，上，右)；外边距 omi
dev.off()

PredictEnvs <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                 "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")
# 循环迭代每个环境因子
for (env_factor in PredictEnvs) {
  pdf_file <- paste0(picture_dir, "/cumulativeplot_", env_factor, ".pdf")
  pdf(file = pdf_file, width = 5, height = 7)
  plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = env_factor, leg.nspecies = 5,
       show.overall = T, legend = T, leg.posn = "topleft", common.scale = T, 
       cex.lab = 1, cex.legend = 0.7, cex.axis = 1, line.ylab = 2,
       par.args = list(mgp = c(2, 0.5, 0), mar = c(0.5, 0.2, 0.2, 1), 
                       omi = c(0.3, 0.7, 0.1, 0.1)))
  dev.off()
}


#显示整体（环境因子bio1-19）组成的累积变化，其中变化发生在梯度上
pdf(file = paste0(picture_dir, "/All_predictorcumulative_common.scale_F.pdf"), 
    width = 8, height = 8)
plot(gf.mod, plot.type = "C", imp.vars = PredictEnvsNames, show.species = F, 
     common.scale = F, cex.axis = 1, cex.lab = 1.2, line.ylab = 0.8, 
     par.args = list(mgp = c(1.7, 0.5, 0), mar = c(3, 2, 0.1, 0.5), 
                     omi = c(0, 0.3, 0.05, 0.1)))
dev.off()


# R2 重要性
pdf(file = paste0(picture_dir, "/All_R2_performance_rank.pdf"), width = 7, height = 7)
plot(gf.mod, plot.type = "Performance", show.names = T, horizontal = T, 
     cex.axis = 1.2, cex.labels = 0.8, line = 2, 
     par.args = list(mgp = c(0, 0.8, 0),
                     mar = c(4, 7, 2, 0.5), omi = c(0, 0, 0.1, 0.1)))
dev.off()

###########################################################################

# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
All_current_xy_envs = fread("extracted_future_data/future_climate_current_O.fragrans.csv")
All_current_xy_envs = as.data.frame(All_current_xy_envs)
# 从 All_current_xy_envs 数据框中筛选出指定列和第一列、第二列，并且只保留这些列中不包含缺失值的行
All_current_xy_envs = All_current_xy_envs[complete.cases(All_current_xy_envs[, c("lon", "lat", PredictEnvs)]), 
                                          c("lon", "lat", PredictEnvs)]
head(All_current_xy_envs)
dim(All_current_xy_envs)

# 将 All_current_xy_envs 数据与梯度森林模型 gf.mod 计算，量化环境梯度
All_grids = cbind(All_current_xy_envs[, c("lon", "lat")], 
                  predict(gf.mod, All_current_xy_envs[, PredictEnvs]))

# 统计包含NA的数量
n <- sum(is.na(All_grids))
n
# 删除包含NA的行，得到没有缺失值的数据
Trns_grid <- na.omit(All_grids)
# 重新统计包含NA的数量
n <- sum(is.na(Trns_grid))
n
# 再次删除包含NA的行，确保数据没有缺失值
Trns_grid <- na.omit(All_grids)
# 重新统计包含NA的数量
n <- sum(is.na(Trns_grid))
n

# 使用主成分分析（PCA）进行降维，选择特定的环境因子进行分析
# 基于梯度森林结果重要值排序和排除自相关筛选出的环境因子
All_PCs <- prcomp(Trns_grid[, PredictEnvs], 
                  center = TRUE, scale. = FALSE)
# 输出PCA分析的摘要信息
summary(All_PCs)

# 为地图绘制设置颜色，使用主成分的分数进行颜色映射
a1 <- All_PCs$x[, 1]
a2 <- All_PCs$x[, 2]
a3 <- All_PCs$x[, 3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
# 归一化颜色值到0-255范围
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
# 创建包含坐标和颜色信息的数据框
grid <- All_current_xy_envs[, c("lon","lat")]
grid$R = r
grid$G = g
grid$B = b
# 获取主成分分析中环境因子的数量
nvs <- dim(All_PCs$rotation)[1]
nvs
vec <- PredictEnvs
lv <- length(vec)
lv
vind <- rownames(All_PCs$rotation) %in% vec
scal <- 60
xrng <- range(All_PCs$x[, 1], All_PCs$rotation[, 1]/scal) * 1.1
yrng <- range(All_PCs$x[, 2], All_PCs$rotation[, 2]/scal) * 1.1

############################  PC散点图  #################################
library(ggplot2)

# 准备主成分得分数据 (PC1 和 PC2)
pc_scores <- as.data.frame(All_PCs$x[, 1:2])
colnames(pc_scores) <- c("PC1", "PC2")

# 添加颜色数据 (RGB)
pc_scores$R <- r
pc_scores$G <- g
pc_scores$B <- b
pc_scores$Color <- rgb(pc_scores$R, pc_scores$G, pc_scores$B, maxColorValue = 255)

# 准备方向向量数据
arrows_data <- as.data.frame(All_PCs$rotation[, 1:2])
colnames(arrows_data) <- c("ArrowX", "ArrowY")
arrows_data$Var <- vec  # 添加环境因子名称
arrows_data[, c("ArrowX", "ArrowY")] <- arrows_data[, c("ArrowX", "ArrowY")] / scal # 对数据框的数值列进行缩放

# 绘制主成分得分散点图
pdf(file = paste0(picture_dir, "/PCplot_PC1_PC2.pdf"), width = 7.5, height = 5.5)
ggplot() +
  # 绘制散点图
  geom_point(data = pc_scores, aes(x = PC1, y = PC2, color = Color), size = 1.2, shape = 16) +
  scale_color_identity() +  # 直接使用提供的颜色
  # 添加方向向量
  geom_segment(data = arrows_data, aes(x = 0, y = 0, xend = ArrowX, yend = ArrowY),
               arrow = arrow(length = unit(0.2, "cm")), color = "black", linewidth = 0.4) +  # 使用 linewidth 替换 size
  # 添加环境因子标签
  geom_text(data = arrows_data, aes(x = ArrowX + jit * sign(ArrowX), 
                                    y = ArrowY + jit * sign(ArrowY), 
                                    label = Var),
            size = 4, color = "black") +
  coord_fixed(ratio = 1) +  # 确保 x 和 y 轴比例一致
  theme_bw() + # 简洁的主题
  # 调整轴标签字体和颜色
  theme(
    axis.title.x = element_text(size = 16, color = "black"),  # x轴标签字体大小和颜色
    axis.title.y = element_text(size = 16, color = "black"), # y轴标签字体大小和颜色
    axis.text = element_text(size = 14, color = "black"),    # x和y轴刻度文字大小
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank()) +  # 去掉次网格线
  labs(x = "Principal Component 1 (PC1)", 
       y = "Principal Component 2 (PC2)", 
       title = NULL) 
dev.off()

############################  PC地图   #################################
library(ggplot2)

# 创建包含经纬度和颜色的绘图数据框
map_data <- as.data.frame(Trns_grid[, c("lon", "lat")])  # 提取经纬度数据
colnames(map_data) <- c("Longitude", "Latitude")
map_data$Color <- rgb(r, g, b, max = 255)  # 将颜色信息加入数据框

# 绘制主成分地图
pdf(file = paste0(picture_dir, "/PCplot_MAP.pdf"), width = 7.5, height = 5.5)
ggplot(map_data, aes(x = Longitude, y = Latitude)) +
  # 绘制点图
  geom_point(aes(color = Color), size = 1e-10) +  # size 控制点的大小
  scale_color_identity() +  # 直接使用 RGB 颜色
  coord_fixed(ratio = 1) +  # 确保经纬度比例一致
  scale_y_continuous(limits = c(19, 35)) +  # 设置 Y 轴范围，
  theme_bw() +  # 使用简洁的主题
  # 设置标题和轴标签
  labs(
    title = NULL,  # 添加主标题
    x = "Longitude", 
    y = "Latitude") +
  # 调整主题样式
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 标题居中加粗
    axis.title = element_text(size = 16, color = "black"),  # 坐标轴 标签字体大小颜色
    axis.text = element_text(size = 14, color = "black"),  # 坐标轴 刻度字体大小颜色
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    panel.background = element_blank())  # 移除背景填充
dev.off()

##########################  保存PC颜色信息 ####################################

# 将RGB颜色转换为ArcGIS使用的格式
greencols=rgb(r,g,b, max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)

# 创建包含颜色信息的数据集
gradients=cbind(Trns_grid[, 1:2],greencols3)
gradients$color=greencols

# 将颜色数据保存为CSV文件
write.csv(gradients,file = paste0(result_dir, "/all_gradients4arcgis.csv"),row.names=F,quote=F)

# 将 All_current_xy_envs 数据与梯度森林模型 gf.mod 计算，量化环境梯度
All_grids = cbind(All_current_xy_envs[, c("lon", "lat")], 
                  predict(gf.mod, All_current_xy_envs[, PredictEnvs]))


##################### 计算未来场景遗传偏移 LOCAL GENETIC OFFSET #####################

# 读取未来气候数据
All_future_xy_envs <- read.csv("extracted_future_data/future_climate_ssp245_2041-2060_O.fragrans.csv")

# 直接使用逻辑向量来筛选数据，并删除包含NA的行
All_future_xy_envs <- All_future_xy_envs[complete.cases(All_future_xy_envs[, PredictEnvs]), 
                                         c("lon", "lat", PredictEnvs)]

head(All_future_xy_envs)
# 使用梯度森林模型计算环境梯度（未来）
All_future_xy_envs_cbind = cbind(All_future_xy_envs[, c("lon","lat")], 
                                 predict(gf.mod, All_future_xy_envs[, PredictEnvs]))

ncol(All_future_xy_envs_cbind)  ## 应当为预测环境因子数量+2 (环境因子数量+经纬度)
# 计算遗传偏移 所有预测环境因子的遗传偏移
genOffsetAll <- sqrt(rowSums((All_future_xy_envs_cbind[, 3:ncol(All_future_xy_envs_cbind)] - All_future_xy_envs_cbind[, 3:ncol(All_future_xy_envs_cbind)])^2))

# 将遗传偏移值合并到坐标数据中
Offset=cbind(All_future_xy_envs_cbind[, c("lon","lat")], genOffsetAll)

# 修改列名为“offset”
colnames(Offset)[3] <-"offset"

if (!dir.exists("Local_Genetic_Offset")) dir.create("Local_Genetic_Offset")

# 保存遗传偏移结果为CSV文件，以便在ArcGIS中使用
write.csv(Offset, file = paste0(Local_GO_dir, "/Local_Genetic_Offset_ssp245_2041_2060.csv"), quote = FALSE, row.names = FALSE)


#################  循环计算多个未来场景 Local Genetic Offset 遗传偏移 #################   
# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
All_current_xy_envs = fread("extracted_future_data/future_climate_current_O.fragrans.csv")
All_current_xy_envs = as.data.frame(All_current_xy_envs)
# 从 All_current_xy_envs 数据框中筛选出指定列和第一列、第二列，并且只保留这些列中不包含缺失值的行
All_current_xy_envs = All_current_xy_envs[complete.cases(All_current_xy_envs[, c("lon", "lat", PredictEnvs)]), 
                                          c("lon", "lat", PredictEnvs)]
head(All_current_xy_envs)
dim(All_current_xy_envs)

# 将 All_current_xy_envs 数据与梯度森林模型 gf.mod 计算，量化环境梯度
All_grids = cbind(All_current_xy_envs[, c("lon", "lat")], 
                  predict(gf.mod, All_current_xy_envs[, PredictEnvs]))

# 定义时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

# 确保存放结果的目录存在
if (!dir.exists("Local_Genetic_Offset")) {
  dir.create("Local_Genetic_Offset")
}

# 循环遍历每个时期
for (period in periods) {
  # 读取未来气候数据
  file_name <- paste0("extracted_future_data/future_climate_", period, "_O.fragrans.csv")
  future_data <- read.csv(file_name)
  # 筛选数据，并删除包含NA的行
  future_data <- future_data[complete.cases(future_data[, PredictEnvs]), c("lon", "lat", PredictEnvs)]
  # 使用梯度森林模型进行环境梯度预测
  future_data_pred <- cbind(future_data[,c("lon", "lat")], 
                            predict(gf.mod, future_data[, PredictEnvs]))
  # 计算遗传偏移
  genOffsetAll <- sqrt(rowSums((future_data_pred[, 3:ncol(future_data_pred)] - 
                                  All_grids[, 3:ncol(All_grids)])^2))
  # 合并遗传偏移值到坐标数据中
  Offset <- cbind(future_data[, c("lon", "lat")], genOffset = genOffsetAll)
  # 将结果写入 CSV 文件
  output_file <- paste0("Local_Genetic_Offset/Local_Genetic_Offset_", period, ".csv")
  write.csv(Offset, file = output_file, quote = FALSE, row.names = FALSE)
}

##############################################  FORWARD & REVERSE GENETIC OFFSET ##########################################
# 加载必需的包
require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)

# 选择预测所用的环境因子
PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")


######################### 计算正向遗传偏移 ForwardOffset  ######################## 
###  FORWARD GENETIC OFFSET 

# 读取未来气候数据
FutureEnvData <- read.csv("extracted_future_data/future_climate_ssp245_2041-2060_O.fragrans.csv")
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n
# 直接使用逻辑向量来筛选数据，并删除包含NA的行
FutureEnvData <- FutureEnvData[complete.cases(FutureEnvData[, c("lon", "lat", PredictEnvs)]), 
                               c("lon", "lat", PredictEnvs)]
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n

# 使用梯度森林模型转换未来气候数据
FutureEnvDataGF <- data.frame(FutureEnvData[, c("lon", "lat")], predict(gf.mod, FutureEnvData[, PredictEnvs]))

# 使用梯度森林模型转换当前气候数据
# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
All_Current_XY_Envs = read.csv("extracted_future_data/future_climate_current_O.fragrans.csv")
All_Current_XY_Envs = All_Current_XY_Envs[complete.cases(All_Current_XY_Envs[, c("lon", "lat", PredictEnvs)]), 
                                          c("lon", "lat", PredictEnvs)]
n <- sum(is.na(All_Current_XY_Envs))
n
dim(All_Current_XY_Envs)

popDatGF <- data.frame(All_Current_XY_Envs, xy=TRUE, na.rm=TRUE)
popDatGF <- data.frame(All_Current_XY_Envs[, c("lon", "lat")], predict(gf.mod, All_Current_XY_Envs[, PredictEnvs]))
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))
length(popDatGF)

# 正向向遗传偏移计算
cl <- makeCluster(60)
registerDoParallel(cl)

forwardOffsetGF <- foreach(i = 1:length(popDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  print(i)
  onePopGF <- popDatGF[[i]]
  combinedDatGF <- FutureEnvDataGF[,c("lon","lat")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,PredictEnvs], FutureEnvDataGF[,PredictEnvs]))
  coordGF <- onePopGF[,c("lon","lat")]
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  minValGF <- minCoordsGF$gfOffset
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("lon","lat")]
  bearGF <- bearing(coordGF, minPtGF)
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]], y2=minPtGF[[2]])
}

stopCluster(cl)
forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
write.csv(forwardOffsetGF, file = paste0(Forward_GO_dir, "/Forward_Genetic_Offset_ssp245_2041_2060.csv"), row.names = FALSE)

######################### 计算正向遗传偏移 ForwardOffset  限制距离 50KM ######################## 
# 读取未来气候数据
FutureEnvData <- read.csv("extracted_future_data/future_climate_ssp245_2061-2080_O.fragrans.csv")
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n
# 直接使用逻辑向量来筛选数据，并删除包含NA的行
FutureEnvData <- FutureEnvData[complete.cases(FutureEnvData[, c("lon", "lat", PredictEnvs)]), 
                               c("lon", "lat", PredictEnvs)]
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n

# 使用梯度森林模型转换未来气候数据
FutureEnvDataGF <- data.frame(FutureEnvData[, c("lon", "lat")], predict(gf.mod, FutureEnvData[, PredictEnvs]))

# 使用梯度森林模型转换当前气候数据
# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
All_Current_XY_Envs = read.csv("extracted_future_data/future_climate_current_O.fragrans.csv")
All_Current_XY_Envs = All_Current_XY_Envs[complete.cases(All_Current_XY_Envs[, c("lon", "lat", PredictEnvs)]), 
                                          c("lon", "lat", PredictEnvs)]
n <- sum(is.na(All_Current_XY_Envs))
n
dim(All_Current_XY_Envs)

popDatGF <- data.frame(All_Current_XY_Envs, xy=TRUE, na.rm=TRUE)
popDatGF <- data.frame(All_Current_XY_Envs[, c("lon", "lat")], predict(gf.mod, All_Current_XY_Envs[, PredictEnvs]))
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))
length(popDatGF)

cl <- makeCluster(60)
registerDoParallel(cl)
forwardOffsetGF <- foreach(i = 1:length(popDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
  onePopGF <- popDatGF[[i]]
  combinedDatGF <- FutureEnvDataGF[,c("lon","lat")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,PredictEnvs], FutureEnvDataGF[,PredictEnvs]))
  coordGF <- onePopGF[,c("lon","lat")]
  combinedDatGF['dists']=distGeo(p1=coordGF, p2=combinedDatGF[,1:2])
  combinedDatGF<-combinedDatGF[combinedDatGF['dists']<50000,]
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dist == min(minCoordsGF$dists)),]
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y ==coordGF$y),"gfOffset"]
  minValGF <- minCoordsGF$gfOffset
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[,c("lon","lat")]
  bearGF <- bearing(coordGF, minPtGF)
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, forwardOffset=minValGF, predDist=toGoGF, bearing=bearGF,x2=minPtGF[[1]],y2=minPtGF[[2]])
}
stopCluster(cl)
write.csv( do.call(rbind, forwardOffsetGF),paste0("./future_ssp245_2061-2080_50km_ForwardOffsetGF.csv"), row.names=FALSE)


####################### 计算反向遗传偏移 ReverseOffset  ######################## 
# 反向遗传偏移计算
# 读取未来气候数据
FutureEnvData <- read.csv("extracted_future_data/future_climate_ssp245_2061-2080_O.fragrans.csv")
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n
# 直接使用逻辑向量来筛选数据，并删除包含NA的行
FutureEnvData <- FutureEnvData[complete.cases(FutureEnvData[, c("lon", "lat", PredictEnvs)]), 
                               c("lon", "lat", PredictEnvs)]
dim(FutureEnvData)
n <- sum(is.na(FutureEnvData))
n

# 使用梯度森林模型转换未来气候数据
FutureEnvDataGF <- data.frame(FutureEnvData[, c("lon", "lat")], predict(gf.mod, FutureEnvData[, PredictEnvs]))

# 使用梯度森林模型转换当前气候数据
# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
All_Current_XY_Envs = read.csv("extracted_future_data/future_climate_current_O.fragrans.csv")
All_Current_XY_Envs = All_Current_XY_Envs[complete.cases(All_Current_XY_Envs[, c("lon", "lat", PredictEnvs)]), 
                                          c("lon", "lat", PredictEnvs)]
n <- sum(is.na(All_Current_XY_Envs))
n
dim(All_Current_XY_Envs)

popDatGF <- data.frame(All_Current_XY_Envs, xy=TRUE, na.rm=TRUE)
popDatGF <- data.frame(All_Current_XY_Envs[, c("lon", "lat")], predict(gf.mod, All_Current_XY_Envs[, PredictEnvs]))
###  popDatGF <- split(popDatGF, seq(nrow(popDatGF)))   # 注意这里反向遗传偏移计算 不需要split
dim(popDatGF)

cl <- makeCluster(60) 
registerDoParallel(cl) 
reverseOffsetGF <- foreach(i=1:nrow(FutureEnvDataGF), .packages=c("fields", "gdm", "geosphere")) %dopar%{
  print(i)
  onePopGF <- FutureEnvDataGF[i,]
  combinedDatGF <- popDatGF[, c("lon","lat")]
  combinedDatGF["gfOffset"] <- c(rdist(onePopGF[, PredictEnvs], popDatGF[, PredictEnvs]))
  coordGF <- onePopGF[, c("lon","lat")]
  minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
  minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[, 1:2])
  minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
  minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF), 1),]
  offsetGF <- combinedDatGF[which(combinedDatGF$x == coordGF$x & combinedDatGF$y == coordGF$y), "gfOffset"]
  minValGF <- minCoordsGF$gfOffset
  toGoGF <- minCoordsGF$dists
  minPtGF <- minCoordsGF[, c("lon","lat")]
  bearGF <- bearing(coordGF, minPtGF)
  outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]], local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]], y2=minPtGF[[2]])
}

stopCluster(cl) # 停止并行集群
reverseOffsetGF <- do.call(rbind, reverseOffsetGF) # 合并结果
write.csv(reverseOffsetGF, paste0("./future_GF_ssp245_2061-2080_ReverseOffsetGF.csv"), row.names=FALSE) # 保存结果为CSV文件


# 结果数据框包含的列：
# x1/y1: 焦点坐标
# local: 局部偏移
# reverseOffset: 逆向偏移
# predDist: 到逆向偏移地点的距离
# bearing: 到逆向偏移地点的方位角
# x2/y2: 逆向偏移地点的坐标
