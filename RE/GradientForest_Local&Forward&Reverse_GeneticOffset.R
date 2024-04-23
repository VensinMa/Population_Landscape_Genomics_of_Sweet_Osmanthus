# 加载所需的R包
library(gradientForest)
library(MaizePal)
library(data.table)
library(dplyr)
library(tidyverse)

getwd()

# 设置导出结果目录和图片目录
result_dir <- "final_result"
picture_dir <- "final_picture"

# 检查输出目录是否存在，如果不存在，则创建
if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(picture_dir)) dir.create(picture_dir)

# 导入群体适应性位点等位基因频率的数据
PopsMaf <- fread(file = "input/1861SNP_32pops_combined_maf.csv")
setnames(PopsMaf, old = "V1", new = "pop") # 第一列是群体名称但没有列名 需添加
setDF(PopsMaf) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(PopsMaf) <- PopsMaf[[1]] # 假设我们想要将第一列作为行名
PopsMaf <- PopsMaf[,-1] # 移除已经设置为行名的列
PopsMaf <- PopsMaf[order(rownames(PopsMaf)), ]
PopsMaf <- PopsMaf %>%
  select(where(~ !any(is.na(.))))
dim(PopsMaf)

## 选择保留的预测环境因子
PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

# 读取当前环境数据
CurrentEnvs <- read.csv('input/32pop_means_env_vars.csv', 
                        header = T, row.names = 1)
CurrentEnvs <- CurrentEnvs[, PredictEnvs]

PredictEnvsNames = colnames(CurrentEnvs)
print(PredictEnvsNames)

# 计算maxLevel，用于树的最大深度
maxLevel <- log2(0.368 * nrow(PopsMaf) / 2)
Envs_Maf = cbind(CurrentEnvs[, PredictEnvsNames], PopsMaf)

#load(file = "gf.mod.RData")

# 构建梯度森林模型，将环境数据和等位基因频率合并在一起
gf.mod <- gradientForest(Envs_Maf,
                         predictor.vars = colnames(CurrentEnvs),# 指定了用作预测变量的列名，即环境数据的列名
                         response.vars = colnames(PopsMaf),# 指定了用作响应变量的列名，即等位基因频率数据的列名
                         ntree = 1000, maxLevel = maxLevel, trace = T,
                         corr.threshold = 0.50,  nbin = 1001, check.names = FALSE)

# 保存gf_all_SNP到一个文件
#save(gf.mod, file = "gf.mod.RData")

############################## 绘图 #########################################

# 绘制特征重要性图，不同颜色表示特征的整体重要性
#生成空的PDF文件 #生成重要值排序 #保存生成的结果
pdf(file = paste0(picture_dir, "/gf.mod.Importance.pdf"), width = 8, height = 8) 
plot(gf.mod, plot.type = "Overall.Importance", 
     col = c(rep("grey",8), MaizePal::maize_pal("HighlandMAGIC", 3)),
     las = 2, cex.names = 0.8) 
dev.off() 

# 将梯度森林模型的各种输出数据写入文本文件 
# 包括响应变量 Y、预测变量 X、重要性值 imp.rsq、模型结果 result等
write.table(gf.mod$Y, file = paste0(result_dir, "/gf.mod.txt"))
write.table(gf.mod$X, file = paste0(result_dir, "/all_gf_X.txt"))
write.table(gf.mod$imp.rsq, file = paste0(result_dir, "/all_gf_impRsq.txt"))
write.table(gf.mod$result, file = paste0(result_dir, "/all_gf_result.txt"))
write.table(gf.mod$res.u, file = paste0(result_dir, "/all_gf_res_u.txt"))
write.table(gf.mod$res, file = paste0(result_dir, "/all_gf_res.txt"))

# splits density plots 分割密度图
pdf(file = paste0(picture_dir, "/All_splitsdensityplots.pdf"), 
    width = 10, height = 7)
plot(gf.mod, plot.type= "S", imp.vars = PredictEnvsNames, leg.posn = "topleft", 
     cex.legend = 1.1, cex.axis = 1, cex.lab = 1.2, line.ylab = 0.5, 
     par.args = list(mgp = c(1.6, 0.5, 0), mar = c(3.1,1.5,0.1,1)))
dev.off()

# 全部环境因子  图例只标注前5个最具响应性的snp
pdf(file = paste0(picture_dir, "/All_env_factor_cumulativeplot_common.scale_F.pdf"), 
    width = 16, height = 18)
plot(gf.mod, plot.type = "Cumulative.Importance", imp.vars = PredictEnvsNames, 
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

# 设置绘图输出到文件，创建空的PDF文件
pdf(file = paste0(picture_dir, "/All_PCplot02.pdf"), width = 8, height = 7)
# 绘制主成分分析的PC1和PC2
plot((All_PCs$x[, 1:2]), xlim = xrng, ylim = yrng, pch = ".", cex = 7, 
     col = rgb(r, g, b, max = 255), asp = 1)
# 添加主成分分析的方向向量
arrows(rep(0, lv), rep(0, lv), All_PCs$rotation[, 1]/scal, All_PCs$rotation[, 2]/scal, 
       length = 0.1)
# 在图中添加环境因子的标签
jit <- 0.0015
text(All_PCs$rotation[, 1]/scal + jit*sign(All_PCs$rotation[, 1]), 
     All_PCs$rotation[, 2]/scal + jit*sign(All_PCs$rotation[, 2]), 
     labels = vec)
# 结束并保存绘图
dev.off()

# 设置绘图输出到文件，创建空的PDF文件
pdf(file = paste0(picture_dir, "All_PCplot_Map2.pdf"), width = 8, height = 7)
# 使用梯度森林模型进行预测
green.pred <- predict(gf.mod, All_current_xy_envs[, PredictEnvs])
# 绘制地图
plot(Trns_grid[, c("lon", "lat")], pch = 15, cex = 1, asp = 1, 
     col = rgb(r, g, b, max = 255), main = "SNP turnover in O. fragrans")
# 结束并保存绘图
dev.off()


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


##################### 计算未来场景genetic_offset 遗传偏移 #####################

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

if (!dir.exists("Genetic_Offset")) dir.create("Genetic_Offset")

# 保存遗传偏移结果为CSV文件，以便在ArcGIS中使用
write.csv(Offset, "Genetic_Offset/ssp245_2041_2060_genetic_offset111111.csv", quote=F, row.names=F)


#################  循环计算多个未来场景genetic_offset 遗传偏移 #################   
# 定义时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

# 确保存放结果的目录存在
if (!dir.exists("Final_Genetic_Offset")) {
  dir.create("Final_Genetic_Offset")
}

# 循环遍历每个时期
for (period in periods) {
  # 读取未来气候数据
  file_name <- paste0("extracted_future_data/future_climate_", period, "_O.fragrans.csv")
  future_data <- read.csv(file_name)
  # 筛选数据，并删除包含NA的行
  future_data <- future_data[complete.cases(future_data[, PredictEnvs]), c("lon", "lat", PredictEnvs)]
  # 使用梯度森林模型进行环境梯度预测
  future_data_pred <- cbind(future_data[,c("lon","lat")], 
                            predict(gf.mod, future_data[, PredictEnvs]))
  # 计算遗传偏移
  genOffsetAll <- sqrt(rowSums((future_data_pred[, 3:ncol(future_data_pred)] - 
                                  All_grids[, 3:ncol(All_grids)])^2))
  # 合并遗传偏移值到坐标数据中
  Offset <- cbind(future_data_pred[,c("lon","lat")], genOffsetAll)
  colnames(Offset)[3] <-"offset"
  # 保存遗传偏移结果为CSV文件
  output_file_name <- paste0("Final_Genetic_Offset/", period, "_genetic_offset_11ev.csv")
  write.csv(Offset, output_file_name, quote=FALSE, row.names=FALSE)
}

############################################## LOCAL FORWARD REVERSE ##########################################
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

# 读取当前群体位置数据
pops <- read.csv("input/32pop_means_env_vars.csv")
pops <- pops[, c(1, 2, 3)]
colnames(pops) <- c('pop', 'lon', 'lat')

# 选择预测所用的环境因子
PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")


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
dim(popDatGF)

####################### 计算正向遗传偏移 ForwardOffset  ######################## 
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
write.csv(forwardOffsetGF, paste0("./future_GF_ssp245_2041-2060_ForwardOffsetGF.csv"), row.names=FALSE)

####################### 计算反向遗传偏移 ReverseOffset  ######################## 
# 反向遗传偏移计算
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
write.csv(reverseOffsetGF, paste0("./future_GF_ssp245_2041-2060_ReverseOffsetGF.csv"), row.names=FALSE) # 保存结果为CSV文件

# 结果数据框包含的列：
# x1/y1: 焦点坐标
# local: 局部偏移
# reverseOffset: 逆向偏移
# predDist: 到逆向偏移地点的距离
# bearing: 到逆向偏移地点的方位角
# x2/y2: 逆向偏移地点的坐标

##########################################################################################################################################################################################################################
# 选择预测所用的环境因子
PredictEnvs = c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")
# 读取未来气候数据
FutureEnvData <- read.csv("extracted_future_data/future_climate_ssp585_2041-2060_O.fragrans.csv")
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
dim(popDatGF)

####################### 计算正向遗传偏移 ForwardOffset  ######################## 
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

write.csv(forwardOffsetGF, paste0("./future_GF_ssp245_2041-2060_ForwardOffsetGF.csv"), row.names=FALSE)

####################### 计算反向遗传偏移 ReverseOffset  ######################## 
# 反向遗传偏移计算
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
write.csv(reverseOffsetGF, paste0("./future_GF_ssp245_2041-2060_ReverseOffsetGF.csv"), row.names=FALSE) # 保存结果为CSV文件

# 结果数据框包含的列：
# x1/y1: 焦点坐标
# local: 局部偏移
# reverseOffset: 逆向偏移
# predDist: 到逆向偏移地点的距离
# bearing: 到逆向偏移地点的方位角
# x2/y2: 逆向偏移地点的坐标
