library(gradientForest)
library(MaizePal)
library(data.table)
library(gdm)
library(dplyr)
library(tidyverse)
library(raster)
library(fields)
library(geosphere)
library(gdm)
library(foreach)
library(parallel)
library(doParallel)
library(gradientForest)
library(fields)
library(ggplot2)
library(sf)

# 导入群体适应性位点等位基因频率数据
PopsMaf <- fread(file = "input/GF_PopsMaf.csv")
setnames(PopsMaf, old = "V1", new = "pop") # 第一列是群体名称但没有列名 将默认列名V1更改为pop
setDF(PopsMaf) # 将data.table转换为data.frame，因为行名是data.frame的特性
rownames(PopsMaf) <- PopsMaf[[1]] # 将第一列作为行名
PopsMaf <- PopsMaf[,-1] # 移除已经设置为行名的列
PopsMaf <- PopsMaf[order(rownames(PopsMaf)), ]
PopsMaf <- PopsMaf %>%
  dplyr::select(where(~ !any(is.na(.))))  # 去除缺失位点
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

gfsnp_final <- gf.mod
DH_current_bio <- read.csv("input/future_climate_current_O.fragrans.csv")
imp.vars <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
              "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")

DH_current_bio_clean <- na.omit(DH_current_bio)

df_pca <- cbind(DH_current_bio_clean[, c("lon", "lat")], predict(gfsnp_final, DH_current_bio_clean[, imp.vars]))

PC <- prcomp(df_pca[,imp.vars])  
pcx <- PC$x

#Assign PC
pc1 <- pcx[,1]  
pc2 <- pcx[,2]
pc3 <- pcx[,3]

#define RGB color palette (your choice)
r <- pc2
g <- pc3+pc1-pc2
b <- pc3-pc2
r <- (r-min(r))/(max(r)-min(r))
g <- (g-min(g))/(max(g)-min(g))
b <- (b-min(b))/(max(b)-min(b))

summary(r)
summary(g)
summary(b)

#Biplot
plot(pcx[,1:2],pch = ".", cex = 1,col = rgb(r,g,b),asp = 1)
vec <- imp.vars #Here we chose some important climate variables
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec
class(vind) 
arrow1 <- PC$rotation[c(1,2,3,4,5,6,7,8),1]  
arrow2 <- PC$rotation[c(1,2,3,4,5,6,7,8),2]
arrow_scale <- 30  #set a scale for the length of the arrow
plot(pcx[,1:2],pch = ".", cex = 2,col = rgb(r,g,b),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
jit <- 0.0012 #distance between the text to the arrow
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)

plot(df_pca[,c("lon","lat")],pch = ".",cex =0.1 ,asp = 1,col = rgb(r,g,b))#Plot spatail map with the GF-predicted results

library("knitr")
library("ggplot2")
library("factoextra")
library("cluster")
dir.create("output")
dim(pcx)  
PC <- prcomp(df_pca[, imp.vars])  # 只选取重要变量列进行 PCA
pcx <- PC$x
dim(pcx)

set.seed(123)
sample_size <- min(30000, nrow(pcx))  
pcx_sampled <- pcx[sample(1:nrow(pcx), size = sample_size), ]
dim(pcx_sampled)  # 检查采样后的数据维度

k.max <- min(10, nrow(pcx) - 1)
f <- fviz_nbclust(pcx_sampled, clara, method = "wss", k.max = k.max)
#f <- fviz_nbclust(pcx_sampled, clara, method = "wss", k.max = 16)
###f <- fviz_nbclust(pcx,clara, method = "wss", k.max = 10)
variation <- f$data$y
write.csv(variation, file = "output/variation_all.csv", row.names = FALSE)

# 定义聚类数
k_values <- 1:length(variation)

# 绘制肘部图
plot(k_values, variation, type = "b", pch = 19, frame = FALSE,
     xlab = "聚类数 (K)", ylab = "变异度 (Variation)",
     main = "肘部法则")
# 计算斜率变化
slopes <- diff(variation)
# 计算斜率的二阶差分
second_derivative <- diff(slopes)

# 找到肘点
elbow_point_index <- which.max(second_derivative) + 1  # 加1因为 diff 函数减少了长度
cat("肘点对应的聚类数为:", elbow_point_index, "\n")

library(cluster)
# 5.1 Cluster points into zones
#-------------------------
ncl <- 2 #number of zones determined based on section 4 results
clPCs <- clara(pcx,ncl,sampsize=10000)

#set up the medoid color palette
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

summary(medcolR)
summary(medcolG)
summary(medcolB)

#PCA biplot into groups
plot(pcx[,1:2], pch = ".", cex = 1, 
     col=rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)

# GF-based seed and breeding zones
#------------------------- 
colnames(df_pca)

plot(df_pca[, c("lon","lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]), 
     main = "", 
     xlim = range(df_pca$lon),        # 修改为数据的范围
     ylim = range(df_pca$lat))         # 修改为数据的范围
# 添加图例
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = rgb(medcolR, medcolG, medcolB))
# 标记聚类中值点
points(df_pca[clPCs$i.med, c("lon","lat")], 
       pch = as.character(seq(1, ncl)))
mix <- rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering])#mix is the colors for each location
summary(mix)
str(mix)
unique(mix)
mix2 <- gsub("#836A69", "1", mix) #extract color code
mix3 <- gsub("#92B970" , "2", mix2)

mix9 <- as.numeric(mix3)
mix9
unique(mix9)
mixcluster <- df_pca[,c("lon","lat")]
mixcluster$color <- mix9
#save in CSV format
write.csv(mixcluster,"output/rzone2.csv")

########################
ncl <- 4  # 设置为 4 个聚类
clPCs <- clara(pcx, ncl, sampsize = 10000)

# 设置颜色
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

# 绘制 PCA 双变量图
plot(pcx[, 1:2], 
     pch = ".", 
     cex = 1, 
     col = rgb(medcolR[clPCs$clustering], 
               medcolG[clPCs$clustering], 
               medcolB[clPCs$clustering]), 
     asp = 1)
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
text(arrow1 / arrow_scale + jit * sign(arrow1), 
     arrow2 / arrow_scale + jit * sign(arrow2), 
     labels = vec)

# 绘制地理分布图
plot(df_pca[, c("lon", "lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = rgb(medcolR[clPCs$clustering], 
               medcolG[clPCs$clustering], 
               medcolB[clPCs$clustering]), 
     main = "", 
     xlim = range(df_pca$lon, na.rm = TRUE), 
     ylim = range(df_pca$lat, na.rm = TRUE))

# 添加图例
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = rgb(medcolR, medcolG, medcolB))

# 标记中值点
points(df_pca[clPCs$i.med, c("lon", "lat")], 
       pch = as.character(seq(1, ncl)))

# 处理颜色代码
mix <- rgb(medcolR[clPCs$clustering],
           medcolG[clPCs$clustering],
           medcolB[clPCs$clustering])
mix2 <- gsub("#953E4E", "1", mix)
mix3 <- gsub("#95A56A", "2", mix2)
mix4 <- gsub("#BA5546", "3", mix3)  # 替换为实际颜色
mix5 <- gsub("#935377", "4", mix4)  # 替换为实际颜色

mix9 <- as.numeric(mix5)

# 保存结果
mixcluster <- df_pca[, c("lon", "lat")]
mixcluster$color <- mix9
write.csv(mixcluster, "output/rzone4.csv")



##############################
# 设置聚类数
ncl <- 4  # 设置为任意聚类数

# 运行 CLARA 聚类
clPCs <- clara(pcx, ncl, sampsize = 10000)

# 设置动态颜色
color_palette <- c("#836A69", "#92B970", "#6A89CC", "#D97C7C")  # , "#A66AB8","#7CD9C3"添加更多颜色可扩展
if (ncl > length(color_palette)) {
  stop("请为更多的聚类数添加颜色到 color_palette 中")
}

medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

# 绘制 PCA 双变量图
plot(pcx[, 1:2], 
     pch = ".", 
     cex = 1, 
     col = rgb(medcolR[clPCs$clustering], 
               medcolG[clPCs$clustering], 
               medcolB[clPCs$clustering]), 
     asp = 1)
arrows(rep(0, lv), rep(0, lv), arrow1 / arrow_scale, arrow2 / arrow_scale, length = 0.0625)
text(arrow1 / arrow_scale + jit * sign(arrow1), 
     arrow2 / arrow_scale + jit * sign(arrow2), 
     labels = vec)

# 绘制地理分布图
plot(df_pca[, c("lon", "lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = rgb(medcolR[clPCs$clustering], 
               medcolG[clPCs$clustering], 
               medcolB[clPCs$clustering]), 
     main = "", 
     xlim = range(df_pca$lon, na.rm = TRUE), 
     ylim = range(df_pca$lat, na.rm = TRUE))

# 添加图例
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = color_palette[1:ncl])

# 标记中值点
points(df_pca[clPCs$i.med, c("lon", "lat")], 
       pch = as.character(seq(1, ncl)))

# 动态处理颜色代码
mix <- rgb(medcolR[clPCs$clustering], 
           medcolG[clPCs$clustering], 
           medcolB[clPCs$clustering])

# 动态替换颜色为聚类编号
for (i in 1:ncl) {
  mix <- gsub(color_palette[i], as.character(i), mix)
}

mix9 <- as.numeric(mix)

# 保存结果
mixcluster <- df_pca[, c("lon", "lat")]
mixcluster$color <- mix9
write.csv(mixcluster, sprintf("output/rzone%d.csv", ncl))  # 动态生成文件名

###############################################
gf <- gfsnp_final
write.csv(df_pca, file="output/df_pca.csv", row.names=FALSE)
vars <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
          "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")
df_70245 <- read.csv("/home/vensin/Rstudio/RStudio/Workspace/gradientForest_2024_1m/extracted_future_data_1m/future_climate_ssp245_2061-2080_O.fragrans.csv")
df_70245 <- na.omit(df_70245)  # 去除含有 NA 的行

df_7085 <- read.csv("/home/vensin/Rstudio/RStudio/Workspace/gradientForest_2024_1m/extracted_future_data_1m/future_climate_ssp585_2061-2080_O.fragrans.csv")
df_7085<- na.omit(df_7085)

library(gradientForest)
##
Trns_grid_70245 <- cbind(df_70245[,c("lon","lat")], predict(gf,df_70245[,vars]))
write.csv(Trns_grid_70245, file="output/Trns_grid_70245.csv", row.names=FALSE)
##70-585
Trns_grid_7085 <- cbind(df_7085[,c("lon","lat")], predict(gf,df_7085[,vars]))
write.csv(Trns_grid_7085, file="output/Trns_grid_7085.csv", row.names=FALSE)
#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_7085[, vars])
future_matrix245 <- as.matrix(Trns_grid_70245[, vars]) 

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
knn_result245 <- get.knnx(data = future_matrix245, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
min_distances245 <- knn_result245$nn.dist
min_indices245 <- knn_result245$nn.index
# Create a dataframe with results
#SSP585
results_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_7085$lon[min_indices],
  Future_Latitude = Trns_grid_7085$lat[min_indices],
  Seed_Zone = mixcluster$color
)

write.csv(results_df, "output/Seed_Zone_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv", row.names = FALSE)

#SSP245
results245_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances245,
  Future_Longitude = Trns_grid_70245$lon[min_indices245],
  Future_Latitude = Trns_grid_70245$lat[min_indices245],
  Seed_Zone = mixcluster$color
)

write.csv(results245_df, "output/Seed_Zone_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv", row.names = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(sf)
library(raster)

#ssp585_70
two_zones_df <- read.csv("output/Seed_Zone_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_70585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_70585,"output/arrow_data_70585.csv")

China_provinces <- st_read("中国_省.shp") #  本地导入  ne_states没有台湾省
arrow_data <- arrow_data_70585 

ggplot() +
  # 绘制中国地图背景
  # geom_sf(data = China_country, fill = "grey90", color = "black") +  
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.01) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#953E4E",
    "2" = "#95A56A",
    "3" = "#BA5546",
    "4" = "#935377"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 7, color = "black", stroke = 1.2) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3.5, fontface = "bold", color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 7, color = "black", stroke = 1.2) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3.5, fontface = "bold", color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", linewidth = 0.8) +
  # 设置坐标范围
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  # 设置轴标签
  xlab("Longitude") + 
  ylab("Latitude") +
  # 自定义主题
  theme_bw(base_size = 11) +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )

#ssp245_70
two_zones_df <- read.csv("output/Seed_Zone_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_70245 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_70245,"output/arrow_data_70245.csv")

# 获取中国省级边界数据
# China_provinces <- ne_states(country = "China", returnclass = "sf")
# taiwan_province <- ne_states(country = "taiwan", returnclass = "sf")
China_provinces <- st_read("中国_省.shp") #  本地导入  ne_states没有台湾省
arrow_data <- arrow_data_70245 

ggplot() +
  # 绘制中国地图背景
  # geom_sf(data = China_country, fill = "grey90", color = "black") +  
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.01) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#953E4E",
    "2" = "#95A56A",
    "3" = "#BA5546",
    "4" = "#935377"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 7, color = "black", stroke = 1.2) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3.5, fontface = "bold", color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 7, color = "black", stroke = 1.2) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3.5, fontface = "bold", color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.2, "cm")), 
               color = "black", linewidth = 0.8) +
  # 设置坐标范围
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  # 设置轴标签
  xlab("Longitude") + 
  ylab("Latitude") +
  # 自定义主题
  theme_bw(base_size = 11) +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )

##screen -R guihua
##[detached from 103686.guihua]

##screen -R guihua
##[detached from 10034.guihua]



#screen -R guihua 

#环境提取
library(raster)
data <- read.csv(file.choose(),header=TRUE)   # 样本经纬度信息
names(data)
LonLatData <- data[,c(2,3)]  # "longitude","latitude"
files <- list.files("C:\\Users\\CY\\Desktop\\桂花\\current",pattern='tif',full.names=TRUE)  # 环境数据所在的路径，tif和asc格式提取的结果差不多，就是小数点后的数字有一点点不一样，不影响后续的分析。
files <- list.files("C:\\Users\\CY\\Desktop\\桂花\\SSP245_70",pattern='asc',full.names=TRUE)  # 环境数据所在的路径，tif和asc格式提取的结果差不多，就是小数点后的数字有一点点不一样，不影响后续的分析。
files <- list.files("C:\\Users\\CY\\Desktop\\桂花\\SSP585_70",pattern='asc',full.names=TRUE)  # 环境数据所在的路径，tif和asc格式提取的结果差不多，就是小数点后的数字有一点点不一样，不影响后续的分析。
files <- list.files("C:\\Users\\CY\\Desktop\\桂花\\SSP245_90",pattern='asc',full.names=TRUE)  # 环境数据所在的路径，tif和asc格式提取的结果差不多，就是小数点后的数字有一点点不一样，不影响后续的分析。
files <- list.files("C:\\Users\\CY\\Desktop\\桂花\\SSP585_90",pattern='asc',full.names=TRUE)  # 环境数据所在的路径，tif和asc格式提取的结果差不多，就是小数点后的数字有一点点不一样，不影响后续的分析。



Grids <- raster::stack(files)
Variablesdata <- raster::extract(Grids,LonLatData)
Outfile <- as.data.frame(cbind("Fucusdistichus", LonLatData,Variablesdata))
colnames(Outfile) <- c("species","longitude","latitude", colnames(Variablesdata))
write.csv(Outfile, file = "   .csv")





################################################################################################################################################
###############重新提取

提取点
# 获取数据集中点的数量
num_points <- nrow(sf_data)

# 确保数据点数量足够
if (num_points < 90000) {
  stop("数据集中点的数量少于 90,000")
}

# 
random_points <- sf_data %>%
  sample_n(90000)

# 提取经纬度信息
# 假设你的坐标系是 WGS 84 (EPSG:4326)
random_points <- st_transform(random_points, crs = 4326)  # 转换为 WGS 84 坐标系

# 提取经纬度并转换为数据框
coords <- st_coordinates(random_points)
random_points_df <- data.frame(longitude = coords[,1], latitude = coords[,2])

# 确认抽样后的数据框行数
n_random_points <- nrow(random_points_df)

# 导出为 CSV 文件
output_file <- "extracted_points.csv"  # 输出文件名
write.csv(random_points_df, output_file, row.names = FALSE)

# 输出成功信息
print(paste("成功提取了", n_random_points, "个点，并已保存为", output_file))

# 加载必要的包
library(dplyr)

# 读取已生成的 CSV 文件
input_file <- "extracted_points.csv"  # 替换为你的 CSV 文件路径
points_data <- read.csv(input_file)

# 确保数据点数量足够
num_points <- nrow(points_data)
if (num_points < 100000) {
  stop("数据集中点的数量少于 100,000")
}

# 随机提取 100,000 个点
random_sample <- points_data %>%
  sample_n(100000)

# 导出为新的 CSV 文件
output_file <- "random_sampled_points.csv"  # 输出文件名
write.csv(random_sample, output_file, row.names = FALSE)

# 输出成功信息
print(paste("成功从", input_file, "中随机提取了 100,000 个点，并已保存为", output_file))







重新开始
gfsnp_final <- gf.mod
DH_current_bio <- read.csv("current2.csv")
imp.vars <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12","BIO15", "BIO17", "BIO18")


df_pca <- cbind(DH_current_bio[, c("long", "lat")], predict(gfsnp_final, DH_current_bio[, imp.vars]))

PC <- prcomp(df_pca[,imp.vars])  
pcx <- PC$x

#Assign PC
pc1 <- pcx[,1]  
pc2 <- pcx[,2]
pc3 <- pcx[,3]

#define RGB color palette (your choice)
r <- pc2
g <- pc3+pc1-pc2
b <- pc3-pc2
r <- (r-min(r))/(max(r)-min(r))
g <- (g-min(g))/(max(g)-min(g))
b <- (b-min(b))/(max(b)-min(b))

summary(r)
summary(g)
summary(b)

#Biplot
plot(pcx[,1:2],pch = ".", cex = 1,col = rgb(r,g,b),asp = 1)
vec <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12","BIO15", "BIO17", "BIO18") #Here we chose some important climate variables
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec
class(vind) 
arrow1 <- PC$rotation[c(1,2,3,4,5,6,7,8),1]  
arrow2 <- PC$rotation[c(1,2,3,4,5,6,7,8),2]
arrow_scale <- 30  #set a scale for the length of the arrow
plot(pcx[,1:2],pch = ".", cex = 2,col = rgb(r,g,b),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
jit <- 0.0012 #distance between the text to the arrow
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)

plot(df_pca[,c("long","lat")],pch = ".",cex =0.1 ,asp = 1,col = rgb(r,g,b))#Plot spatail map with the GF-predicted results

library("knitr")
library("ggplot2")
library("factoextra")
library("cluster")
f <- fviz_nbclust(pcx,clara, method = "wss", k.max = 16)###f <- fviz_nbclust(pcx,clara, method = "wss", k.max = 10)
variation <- f$data$y
write.csv(variation, file = "variation_2.csv", row.names = FALSE)

# 定义聚类数
k_values <- 1:length(variation)

# 绘制肘部图
plot(k_values, variation, type = "b", pch = 19, frame = FALSE,
     xlab = "聚类数 (K)", ylab = "变异度 (Variation)",
     main = "肘部法则")
# 计算斜率变化
slopes <- diff(variation)
# 计算斜率的二阶差分
second_derivative <- diff(slopes)

# 找到肘点
elbow_point_index <- which.max(second_derivative) + 1  # 加1因为 diff 函数减少了长度
cat("肘点对应的聚类数为:", elbow_point_index, "\n")
library(cluster)
# 5.1 Cluster points into zones
#-------------------------
ncl <- 2 #number of zones determined based on section 4 results
clPCs <- clara(pcx,ncl,sampsize=10000)

#set up the medoid color palette
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

summary(medcolR)
summary(medcolG)
summary(medcolB)

#PCA biplot into groups
plot(pcx[,1:2], pch = ".", cex = 1, 
     col=rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)

# GF-based seed and breeding zones
#------------------------- 
plot(df_pca[, c("long","lat")], 
     pch = ".", 
     cex = 0.1, 
     asp = 1, 
     col = rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]), 
     main = "", 
     xlim = range(df_pca$long),        # 修改为数据的范围
     ylim = range(df_pca$lat))         # 修改为数据的范围
# 添加图例
legend("bottomleft", 
       as.character(seq(1, ncl)), 
       pch = 15, 
       cex = 1, 
       col = rgb(medcolR, medcolG, medcolB))
# 标记聚类中值点
points(df_pca[clPCs$i.med, c("long","lat")], 
       pch = as.character(seq(1, ncl)))
mix <- rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering])#mix is the colors for each location
summary(mix)
str(mix)
unique(mix)
mix2 <- gsub("#6C8270", "1", mix) #extract color code
mix3 <- gsub("#713769" , "2", mix2)

mix9 <- as.numeric(mix3)
mix9
unique(mix9)
mixcluster <- df_pca[,c("long","lat")]
mixcluster$color <- mix9
#save in CSV format
Write.csv(mixcluster,"rzone2.csv")


移动
gf <- gfsnp_final
write.csv(df_pca, file="df_pca.csv", row.names=FALSE)
vars <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12","BIO15", "BIO17", "BIO18") 
df_70245 <- read.csv("702452.csv")
df_70245 <- df_70245[,2:11] 
names(df_70245)[1] <- "Longitude"
names(df_70245)[2] <- "Latitude"
df_7085 <- read.csv("705852.csv")
df_7085 <- df_7085[,2:11]
names(df_7085)[1] <- "Longitude"
names(df_7085)[2] <- "Latitude"
library(gradientForest)
##
Trns_grid_70245 <- cbind(df_70245[,c("Longitude","Latitude")], predict(gf,df_70245[,vars]))
write.csv(Trns_grid_70245, file="Trns_grid_70245.csv", row.names=FALSE)
##70-585
Trns_grid_7085 <- cbind(df_7085[,c("Longitude","Latitude")], predict(gf,df_7085[,vars]))
write.csv(Trns_grid_7085, file="Trns_grid_7085.csv", row.names=FALSE)
#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_7085[, vars])
future_matrix245 <- as.matrix(Trns_grid_70245[, vars]) 

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix245, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
# Create a dataframe with results
#SSP585
results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_7085$Longitude[min_indices],
  Future_Latitude = Trns_grid_7085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv", row.names = FALSE)

#
knn_result <- get.knnx(data = future_matrix245, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results245_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_70245$Longitude[min_indices],
  Future_Latitude = Trns_grid_70245$Latitude[min_indices]
)

write.csv(results245_df, "Seed_Zone_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv", row.names = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(sf)
library(raster)

#ssp585_70
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_70585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_70585,"arrow_data_70585.csv")
# 获取中国的地理数据
china <- ne_countries(country = "china", returnclass = "sf")

# 获取中国省级边界数据
china_provinces <- ne_states(country = "china", returnclass = "sf")

arrow_data <- arrow_data_70585 
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # 绘制中国地图背景
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # 绘制省级边界
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex =0.01) +
  scale_color_manual(values = c("1" = "#6C8270", "2" = "#713769")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#6C8270", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#713769", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch=c("1","2"), size= 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch=c("1"), size= 5, col = "#6C8270") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch=c("2"), size= 5, col = "#713769") +
  geom_segment(data = arrow_data, aes(x = x+0.53, y = y+0.23, xend = (xend-0.5), yend = (yend-0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size =0.8) +
  coord_sf(xlim = c(90, 125), ylim = c(20, 35), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=14),panel.grid.major = element_blank(), legend.position = "none")


#ssp245_70
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_70245 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_70245,"arrow_data_70245.csv")
# 获取中国的地理数据
china <- ne_countries(country = "china", returnclass = "sf")

# 获取中国省级边界数据
china_provinces <- ne_states(country = "china", returnclass = "sf")

arrow_data <- arrow_data_70245 
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # 绘制中国地图背景
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # 绘制省级边界
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex =0.01) +
  scale_color_manual(values = c("1" = "#6C8270", "2" = "#713769")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#6C8270", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#713769", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch=c("1","2"), size= 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch=c("1"), size= 5, col = "#6C8270") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch=c("2"), size= 5, col = "#713769") +
  geom_segment(data = arrow_data, aes(x = x+0.53, y = y+0.23, xend = (xend-0.5), yend = (yend-0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size =0.8) +
  coord_sf(xlim = c(90, 125), ylim = c(20, 35), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=14),panel.grid.major = element_blank(), legend.position = "none")

##################################90
df_90245 <- read.csv("902452.csv")
df_90245 <- df_90245[,2:11] 
names(df_90245)[1] <- "Longitude"
names(df_90245)[2] <- "Latitude"
df_9085 <- read.csv("905852.csv")
df_9085 <- df_9085[,2:11]
names(df_9085)[1] <- "Longitude"
names(df_9085)[2] <- "Latitude"
library(gradientForest)
##
Trns_grid_90245 <- cbind(df_90245[,c("Longitude","Latitude")], predict(gf,df_90245[,vars]))
write.csv(Trns_grid_90245, file="Trns_grid_90245.csv", row.names=FALSE)
##
Trns_grid_9085 <- cbind(df_9085[,c("Longitude","Latitude")], predict(gf,df_9085[,vars]))
write.csv(Trns_grid_9085, file="Trns_grid_9085.csv", row.names=FALSE)
#install.packages("FNN")
library(FNN)
# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix <- as.matrix(Trns_grid_9085[, vars])
future_matrix245 <- as.matrix(Trns_grid_90245[, vars]) 

# Finding the nearest neighbor
knn_result <- get.knnx(data = future_matrix, query = pca_matrix, k = 1)
#knn_result <- get.knnx(data = future_matrix245, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
# Create a dataframe with results
#SSP585
results_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_7085$Longitude[min_indices],
  Future_Latitude = Trns_grid_7085$Latitude[min_indices]
)

write.csv(results_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv", row.names = FALSE)

#
knn_result <- get.knnx(data = future_matrix245, query = pca_matrix, k = 1)
# Extract minimum distances and their indices
min_distances <- knn_result$nn.dist
min_indices <- knn_result$nn.index
results245_df <- data.frame(
  Original_Longitude = df_pca$long,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances,
  Future_Longitude = Trns_grid_70245$Longitude[min_indices],
  Future_Latitude = Trns_grid_70245$Latitude[min_indices]
)

write.csv(results245_df, "Seed_Zone_Shifts_under_climate_change_2090_ssp245_KNN_method_extended.csv", row.names = FALSE)

library(readxl)
library(dplyr)
library(ggplot2)
library(geosphere)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(sf)
library(raster)

#ssp585_90
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_90585 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_90585,"arrow_data_90585.csv")
# 获取中国的地理数据
china <- ne_countries(country = "china", returnclass = "sf")

# 获取中国省级边界数据
china_provinces <- ne_states(country = "china", returnclass = "sf")

arrow_data <- arrow_data_90585 
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # 绘制中国地图背景
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # 绘制省级边界
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex =0.01) +
  scale_color_manual(values = c("1" = "#6C8270", "2" = "#713769")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#6C8270", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#713769", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch=c("1","2"), size= 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch=c("1"), size= 5, col = "#6C8270") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch=c("2"), size= 5, col = "#713769") +
  geom_segment(data = arrow_data, aes(x = x+0.53, y = y+0.23, xend = (xend-0.5), yend = (yend-0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size =0.8) +
  coord_sf(xlim = c(90, 125), ylim = c(20, 35), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=14),panel.grid.major = element_blank(), legend.position = "none")



#ssp245_90
two_zones_df <- read.csv("Seed_Zone_Shifts_under_climate_change_2090_ssp245_KNN_method_extended.csv")
centroids <- two_zones_df %>%
  group_by(Seed_Zone) %>%
  summarise(
    Orig_Centroid_Lon = mean(Original_Longitude),
    Orig_Centroid_Lat = mean(Original_Latitude),
    Future_Centroid_Lon = mean(Future_Longitude),
    Future_Centroid_Lat = mean(Future_Latitude),
    mean_offset = mean(Min_Distance)
  )
centroids$Distance_km <- mapply(function(lon1, lat1, lon2, lat2) {
  distHaversine(c(lon1, lat1), c(lon2, lat2)) / 1000  # convert from meters to kilometers
}, centroids$Orig_Centroid_Lon, centroids$Orig_Centroid_Lat, centroids$Future_Centroid_Lon, centroids$Future_Centroid_Lat)


arrow_data_90245 <- centroids %>%
  transmute(
    x = Orig_Centroid_Lon, y = Orig_Centroid_Lat,
    xend = Future_Centroid_Lon, yend = Future_Centroid_Lat,
    Seed_Zone)

write.csv(arrow_data_90245,"arrow_data_90245.csv")
# 获取中国的地理数据
china <- ne_countries(country = "china", returnclass = "sf")

# 获取中国省级边界数据
china_provinces <- ne_states(country = "china", returnclass = "sf")

arrow_data <- arrow_data_90245 
ggplot() +
  geom_sf(data = china, fill = "grey90", color = "black") +  # 绘制中国地图背景
  geom_sf(data = china_provinces, fill = NA, color = "black", lwd = 0.5) +  # 绘制省级边界
  geom_point(aes(x = two_zones_df$Future_Longitude, y = two_zones_df$Future_Latitude, col = as.factor(two_zones_df$Seed_Zone)), pch = ".", cex =0.01) +
  scale_color_manual(values = c("1" = "#6C8270", "2" = "#713769")) +
  geom_point(aes(x = centroids[1,]$Orig_Centroid_Lon, y = centroids[1,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#6C8270", color = "black", stroke = 1.2) +
  geom_point(aes(x = centroids[2,]$Orig_Centroid_Lon, y = centroids[2,]$Orig_Centroid_Lat), 
             pch = 21, size = 7, fill = "#713769", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), 
             pch = 21, size = 7, fill = "white", color = "black", stroke=1.2) +
  geom_point(aes(x = centroids$Orig_Centroid_Lon, y = centroids$Orig_Centroid_Lat), pch=c("1","2"), size= 5, col = "white") +
  geom_point(aes(x = centroids[1,]$Future_Centroid_Lon, y = centroids[1,]$Future_Centroid_Lat), pch=c("1"), size= 5, col = "#6C8270") +
  geom_point(aes(x = centroids[2,]$Future_Centroid_Lon, y = centroids[2,]$Future_Centroid_Lat), pch=c("2"), size= 5, col = "#713769") +
  geom_segment(data = arrow_data, aes(x = x+0.53, y = y+0.23, xend = (xend-0.5), yend = (yend-0.5)), 
               arrow = arrow(length = unit(0.2, "cm")), color = "black", size =0.8) +
  coord_sf(xlim = c(90, 125), ylim = c(20, 35), expand = FALSE) +
  xlab("\nLongitude") + ylab("Latitude\n") +
  theme_bw(base_size = 11) +
  theme(plot.background = element_blank(), axis.text = element_text(size=12), axis.title = element_text(size=14),panel.grid.major = element_blank(), legend.position = "none")
