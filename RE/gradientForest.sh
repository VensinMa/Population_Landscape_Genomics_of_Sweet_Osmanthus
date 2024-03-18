# 加载所需的R包
library(gradientForest)
library(MaizePal)
library(data.table)
library(dplyr)

# 设置工作目录
setwd("/public1/guop/mawx/workspace/186sample/gradientForest")

# 检查输出目录是否存在，如果不存在，则创建
if (!dir.exists("result")) dir.create("result")
if (!dir.exists("picture")) dir.create("picture")

# 从 gfData 中提取包含"candSNPs"的列，这些列包含等位基因频率的数据
all_SNP <- fread("/public1/guop/mawx/workspace/186sample/gradientForest/32pop_maf.csv")
setnames(all_SNP, old = "V1", new = "pop") # 第一列是群体名称但没有列名 需添加
# 将第一列设置为行名，并从数据中移除该列
setDF(all_SNP) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(all_SNP) <- all_SNP[[1]] # 假设我们想要将第一列作为行名
all_SNP <- all_SNP[,-1] # 移除已经设置为行名的列
all_SNP <- all_SNP[order(rownames(all_SNP)), ]
all_SNP <- all_SNP %>%
  select(where(~ !any(is.na(.))))

# 环境数据
presClim <- fread('/public1/guop/mawx/workspace/186sample/gradientForest/32pop_means_env_vars.csv')
# 将第一列设置为行名，并从数据中移除该列
setDF(presClim) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(presClim) <- presClim[[1]] # 将第一列作为行名
presClim <- presClim[,-1] # 移除已经设置为行名的列
presClim <- presClim[order(rownames(presClim)), ]
presClim <- presClim[, 3:ncol(presClim)]

bioclimatic = colnames(presClim)
print(bioclimatic)
# 计算maxLevel，用于树的最大深度
maxLevel <- log2(0.368 * nrow(all_SNP) / 2)
cbind_data = cbind(presClim[, bioclimatic], all_SNP)

# 构建梯度森林模型，将环境数据和等位基因频率合并在一起
gf_all_SNP <- gradientForest(cbind_data,
                             predictor.vars = colnames(presClim),# 指定了用作预测变量的列名，即环境数据的列名
                             response.vars = colnames(all_SNP),# 指定了用作响应变量的列名，即等位基因频率数据的列名
                             ntree = 500, maxLevel = maxLevel, trace = T,
                             corr.threshold = 0.50,  nbin =1001, check.names = FALSE)
# 保存 gf_all_SNP 对象到当前工作目录
save(gf_all_SNP, file = "gf_all_SNP.RData")

# 绘制特征重要性图，不同颜色表示特征的整体重要性
#生成空的PDF文件 #生成重要值排序 #保存生成的结果
pdf(file="picture/gf_all_SNP01.pdf") 
plot(gf_all_SNP, plot.type = "Overall.Importance", 
     col=c(rep("grey",15),MaizePal::maize_pal("HighlandMAGIC", 4) ),
     las=2,cex.names=0.8) 
dev.off() 

# 将梯度森林模型的各种输出数据写入文本文件 
# 包括响应变量 Y、预测变量 X、重要性值 imp.rsq、模型结果 result、独特的行 res.u 和合并的数据框 res
write.table(gf_all_SNP$Y, file="result/gf_all_SNP.txt")
write.table(gf_all_SNP$X, file="result/all_gf_X.txt")
write.table(gf_all_SNP$imp.rsq, file="result/all_gf_impRsq.txt")
write.table(gf_all_SNP$result, file="result/all_gf_result.txt")
write.table(gf_all_SNP$res.u, file="result/all_gf_res_u.txt")
write.table(gf_all_SNP$res, file="result/all_gf_res.txt")

#splits density plots 分割密度图
pdf(file="picture/all_splitsdensityplots.pdf")
plot(gf_all_SNP, plot.type="S", imp.vars= bioclimatic, leg.posn="topright", 
     cex.legend=0.4, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, 
     par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))
dev.off()

#图例显示了每个预测因子的前5个最具响应性的snp
pdf(file="picture/all_speciescumulativeplot01.pdf")
plot(gf_all_SNP, plot.type="Cumulative.Importance", imp.vars= bioclimatic, show.overall=T, 
     legend=T,common.scale=T,leg.posn="topleft", leg.nspecies=5, cex.lab=0.7, cex.legend=0.4, 
     cex.axis=0.6, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1),omi=c(0,0.3,0,0)))
dev.off()

#显示整体（环境因子bio1-19）组成的累积变化，其中变化发生在梯度上
pdf(file="picture/all_predictorcumulative.pdf")
plot(gf_all_SNP, plot.type="C", imp.vars= bioclimatic, show.species=F, 
     common.scale=T, cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, par.args=list(mgp=c(1.5, 0.5, 0), 
                                                                             mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))
dev.off()

#R2
pdf(file="picture/all_R2.pdf")
plot(gf_all_SNP, plot.type="P", show.names=F, horizontal=F, cex.axis=1, cex.labels=0.7, line=2.5)
dev.off()

# 从 CSV 文件中读取数据，该文件包含了研究区域的坐标点和气候数据
greengrid = read.csv("wc2.1_30s_bioc_present_xyyjqy.csv")
greengrid = read.csv("extracted_environment_data_香杨研究区域.csv")
# 创建一个逻辑向量，标记包含NA（缺失值）的行
na_rows <- complete.cases(greengrid[, 3:21])
# 使用逻辑向量来筛选数据，删除包含NA的行
greengrid_cleaned <- greengrid[na_rows, 1:21]
# 将清理后的数据重新赋值给greengrid
greengrid = greengrid_cleaned

# 输出 greengrid 的前几行数据
head(greengrid)
# 获取 greengrid 的维度信息，即行数和列数
dim(greengrid)

# 将 greengrid 数据与梯度森林模型 gf_all_SNP 计算，量化环境梯度
all_tgrid = cbind(greengrid[,c("lon","lat")], predict(gf_all_SNP, greengrid[,3:21]))

# 统计包含NA的数量
n <- sum(is.na(all_tgrid))
# 删除包含NA的行，得到没有缺失值的数据
Trns_grid <- na.omit(all_tgrid)
# 重新统计包含NA的数量
n <- sum(is.na(Trns_grid))
# 再次删除包含NA的行，确保数据没有缺失值
Trns_grid <- na.omit(all_tgrid)
# 重新统计包含NA的数量
n <- sum(is.na(Trns_grid))

# 使用主成分分析（PCA）进行降维，选择特定的环境因子进行分析
# 基于梯度森林结果重要值排序和排除自相关筛选出的环境因子
all_PCs <- prcomp(Trns_grid[,c("bio_9","bio_4","bio_13","bio_19","bio_18","bio_14","bio_17")], center=TRUE, scale.=FALSE)
# 输出PCA分析的摘要信息
summary(all_PCs)

# 为地图绘制设置颜色，使用主成分的分数进行颜色映射
a1 <- all_PCs$x[,1]
a2 <- all_PCs$x[,2]
a3 <- all_PCs$x[,3]
r <- a1 + a2
g <- -a2
b <- a3 + a2 - a1
# 归一化颜色值到0-255范围
r <- (r - min(r)) / (max(r) - min(r)) * 255
g <- (g - min(g)) / (max(g) - min(g)) * 255
b <- (b - min(b)) / (max(b) - min(b)) * 255
# 创建包含坐标和颜色信息的数据框
grid <- greengrid[,c("lon","lat")]
grid$R = r
grid$G = g
grid$B = b
# 获取主成分分析中环境因子的数量
nvs <- dim(all_PCs$rotation)[1]
vec <- c("bio_9","bio_4","bio_13","bio_19","bio_18","bio_14","bio_17")
lv <- length(vec)
vind <- rownames(all_PCs$rotation) %in% vec
scal <- 60
xrng <- range(all_PCs$x[,1], all_PCs$rotation[,1]/scal) * 1.1
yrng <- range(all_PCs$x[,2], all_PCs$rotation[,2]/scal) * 1.1

# 设置绘图输出到文件，创建空的PDF文件
pdf(file="picture/all_PCplot02.pdf")
# 绘制主成分分析的PC1和PC2
plot((all_PCs$x[,1:2]), xlim=xrng, ylim=yrng, pch=".", cex=7, col=rgb(r,g,b, max = 255), asp=1)
# 添加主成分分析的方向向量
arrows(rep(0,lv), rep(0,lv), all_PCs$rotation[,1]/scal, all_PCs$rotation[,2]/scal, length = 0.1)
# 在图中添加环境因子的标签
jit <- 0.0015
text(all_PCs$rotation[,1]/scal+jit*sign(all_PCs$rotation[,1]), all_PCs$rotation[,2]/scal+jit*sign(all_PCs$rotation[,2]), labels = vec)
# 结束并保存绘图
dev.off()

# 设置绘图输出到文件，创建空的PDF文件
pdf("picture/Map2.pdf")
# 使用梯度森林模型进行预测
green.pred <- predict(gf_all_SNP, greengrid[,c("bio_9","bio_4","bio_13","bio_19","bio_18","bio_14","bio_17")])
# 绘制地图
plot(Trns_grid[, c("lon", "lat")], pch=15,cex=1.0,asp=1,col=rgb(r,g,b, max=255),main="SNP turnover in Q rugosa")
# 结束并保存绘图
dev.off()


# 将RGB颜色转换为ArcGIS使用的格式
greencols=rgb(r,g,b,max=255)
greencols2=col2rgb(greencols)
greencols3=t(greencols2)

# 创建包含颜色信息的数据集
gradients=cbind(Trns_grid[,1:2],greencols3)
gradients$color=greencols

# 将颜色数据保存为CSV文件
write.csv(gradients,file="result/all_gradients4arcgis03.csv",row.names=F,quote=F)

# 使用梯度森林模型计算环境梯度（当前）
all_tgrid=cbind(greengrid[,c("lon","lat")], predict(gf_all_SNP,greengrid[,3:21]))

##### genetic_offset  遗传偏移

# 读取未来气候数据
future = read.csv("extracted_environment_data_fut_香杨研究区域.csv")

# 创建一个逻辑向量，标记包含NA的行
na_rows_fut <- complete.cases(future[, 3:21])

# 使用逻辑向量来筛选数据
future_cleaned <- future[na_rows, 1:21]

# 现在 future_cleaned 包含了没有NA的数据
future = future_cleaned

# 使用梯度森林模型计算环境梯度（未来）
future_all=cbind(future[,c("lon","lat")], predict(gf_all_SNP,future[,3:21]))

# 计算遗传偏移 所有19个气候因子的遗传偏移
genOffsetAll<-sqrt((future_all[,3]-all_tgrid[,3])^2+
                     (future_all[,4]-all_tgrid[,4])^2+
                     (future_all[,5]-all_tgrid[,5])^2+
                     (future_all[,6]-all_tgrid[,6])^2+
                     (future_all[,7]-all_tgrid[,7])^2+
                     (future_all[,8]-all_tgrid[,8])^2+
                     (future_all[,9]-all_tgrid[,9])^2+
                     (future_all[,10]-all_tgrid[,10])^2+
                     (future_all[,11]-all_tgrid[,11])^2+
                     (future_all[,12]-all_tgrid[,12])^2+
                     (future_all[,13]-all_tgrid[,13])^2+
                     (future_all[,14]-all_tgrid[,14])^2+
                     (future_all[,15]-all_tgrid[,15])^2+
                     (future_all[,16]-all_tgrid[,16])^2+
                     (future_all[,17]-all_tgrid[,17])^2+
                     (future_all[,18]-all_tgrid[,18])^2+
                     (future_all[,19]-all_tgrid[,19])^2+
                     (future_all[,20]-all_tgrid[,20])^2+
                     (future_all[,21]-all_tgrid[,21])^2)

# 将遗传偏移值合并到坐标数据中
Offset=cbind(future_all[,c("lon","lat")],genOffsetAll)

# 修改列名为“offset”
colnames(Offset)[3] <-"offset"

# 保存遗传偏移结果为CSV文件，以便在ArcGIS中使用
write.csv(Offset, "data_2100_24505.csv", quote=F, row.names=T)

# 读取遗传偏移结果
Offset_result = read.csv("data_2100_24505.csv")



