library(gradientForest)
library(MaizePal)
library(data.table)
library(gdm)
library(dplyr)
library(tidyverse)
library(raster)
library(fields)
library(geosphere)
library(foreach)
library(parallel)
library(doParallel)
library(ggplot2)
library(sf)
library(cluster)
library(readxl)
library(rnaturalearth)
library(mapdata)
library(purrr)
library(knitr)
library(factoextra)
library(FNN)

getwd()
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
load(file = "gf.mod.9674.RData")

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

# PCA plot
plot(pcx[,1:2],pch = ".", cex = 1,col = rgb(r,g,b),asp = 1)
vec <- imp.vars #Here we chose some important climate variables
lv <- length(vec)
vind <- rownames(PC$rotation) %in% vec
class(vind) 
arrow1 <- PC$rotation[c(1,2,3,4,5,6,7,8),1]  
arrow2 <- PC$rotation[c(1,2,3,4,5,6,7,8),2]
arrow_scale <- 30  #set a scale for the length of the arrow
plot(pcx[,1:2],pch = ".", cex = 1,col = rgb(r,g,b),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
jit <- 0.0012 #distance between the text to the arrow
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec,cex = 0.5)

plot(df_pca[,c("lon","lat")],pch = ".",cex =0.1 ,asp = 1,col = rgb(r,g,b))#Plot spatail map with the GF-predicted results

dir.create("output")
dim(pcx)  
set.seed(123) 
sample_rows <- sample(1:nrow(df_pca), 40000, replace = FALSE) 
PC <- prcomp(df_pca[sample_rows, imp.vars]) 
pcx_sampled <- PC$x
dim(pcx_sampled)

f <- fviz_nbclust(pcx_sampled, clara, method = "wss", k.max = 16)
#f <- fviz_nbclust(pcx_sampled, clara, method = "wss", k.max = 16)
###f <- fviz_nbclust(pcx,clara, method = "wss", k.max = 10)
variation <- f$data$y
variation
write.csv(variation, file = "output/variation_all.csv", row.names = FALSE)

# 定义聚类数
k_values <- 1:length(variation)

# 绘制肘部图
plot(k_values, variation, type = "b", pch = 19, frame = FALSE,
     xlab = "Number of Clusters (K)", ylab = "Variation",
     main = "Elbow Method")
# 计算斜率变化
slopes <- diff(variation)
# 计算斜率的二阶差分
second_derivative <- diff(slopes)

# 找到肘点
elbow_point_index <- which.max(second_derivative) + 1  # 加1因为 diff 函数减少了长度
cat("肘点对应的聚类数为:", elbow_point_index, "\n")


# 5.1 Cluster points into zones
#-------------------------
ncl <- 2 #number of zones determined based on section 4 results
clPCs <- clara(pcx, ncl, sampsize = 10000)
 
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
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec,cex = 0.6)

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
unique(mix) #[1] "#925372" "#AB8650"
mix2 <- gsub("#925372", "1", mix) #extract color code
mix3 <- gsub("#AB8650" , "2", mix2)

mix9 <- as.numeric(mix3)
mix9
unique(mix9)
mixcluster <- df_pca[,c("lon","lat")]
mixcluster$color <- mix9
#save in CSV format
write.csv(mixcluster,"output/rZone2.csv")

###############################################
gf <- gfsnp_final
write.csv(df_pca, file="output/df_pca.csv", row.names=FALSE)
vars <- c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
          "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")
df_70245 <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024/extracted_future_data/future_climate_ssp245_2061-2080_O.fragrans.csv")
df_70245 <- na.omit(df_70245)  # 去除含有 NA 的行
df_70585 <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024/extracted_future_data/future_climate_ssp585_2061-2080_O.fragrans.csv")
df_70585 <- na.omit(df_70585)
df_90245 <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024/extracted_future_data/future_climate_ssp245_2081-2100_O.fragrans.csv")
df_90245 <- na.omit(df_90245)  # 去除含有 NA 的行
df_90585 <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024/extracted_future_data/future_climate_ssp585_2081-2100_O.fragrans.csv")
df_90585 <- na.omit(df_90585)

##70-245
Trns_grid_70245 <- cbind(df_70245[,c("lon","lat")], predict(gf,df_70245[,vars]))
write.csv(Trns_grid_70245, file="output/Trns_grid_70245.csv", row.names=FALSE)
##70-585
Trns_grid_70585 <- cbind(df_70585[,c("lon","lat")], predict(gf,df_70585[,vars]))
write.csv(Trns_grid_70585, file="output/Trns_grid_70585.csv", row.names=FALSE)
##90-245
Trns_grid_90245 <- cbind(df_90245[,c("lon","lat")], predict(gf,df_90245[,vars]))
write.csv(Trns_grid_90245, file="output/Trns_grid_90245.csv", row.names=FALSE)
##90-585
Trns_grid_90585 <- cbind(df_90585[,c("lon","lat")], predict(gf,df_90585[,vars]))
write.csv(Trns_grid_90585, file="output/Trns_grid_90585.csv", row.names=FALSE)
#install.packages("FNN")

# Convert data frames to matrices for FNN
pca_matrix <- as.matrix(df_pca[, vars])
future_matrix70245 <- as.matrix(Trns_grid_70245[, vars]) 
future_matrix70585 <- as.matrix(Trns_grid_70585[, vars])
future_matrix90245 <- as.matrix(Trns_grid_90245[, vars]) 
future_matrix90585 <- as.matrix(Trns_grid_90585[, vars])

# Finding the nearest neighbor
knn_result70245 <- get.knnx(data = future_matrix70245, query = pca_matrix, k = 1)
knn_result70585 <- get.knnx(data = future_matrix70585, query = pca_matrix, k = 1)
knn_result90245 <- get.knnx(data = future_matrix90245, query = pca_matrix, k = 1)
knn_result90585 <- get.knnx(data = future_matrix90585, query = pca_matrix, k = 1)

# Extract minimum distances and their indices
min_distances70245 <- knn_result70245$nn.dist
min_indices70245 <- knn_result70245$nn.index
min_distances70585 <- knn_result70585$nn.dist
min_indices70585 <- knn_result70585$nn.index

min_distances90245 <- knn_result90245$nn.dist
min_indices90245 <- knn_result90245$nn.index
min_distances90585 <- knn_result90585$nn.dist
min_indices90585 <- knn_result90585$nn.index

# Create a dataframe with results
#SSP70-245
results70245_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances70245,
  Future_Longitude = Trns_grid_70245$lon[min_indices70245],
  Future_Latitude = Trns_grid_70245$lat[min_indices70245],
  Seed_Zone = mixcluster$color
)
write.csv(results70245_df, "output/Seed_Zone2_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv", row.names = FALSE)

#SSP70-585
results70585_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances70585,
  Future_Longitude = Trns_grid_70585$lon[min_indices70585],
  Future_Latitude = Trns_grid_70585$lat[min_indices70585],
  Seed_Zone = mixcluster$color
)
write.csv(results70585_df, "output/Seed_Zone2_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv", row.names = FALSE)

#SSP90-245
results90245_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances90245,
  Future_Longitude = Trns_grid_90245$lon[min_indices90245],
  Future_Latitude = Trns_grid_90245$lat[min_indices90245],
  Seed_Zone = mixcluster$color
)
write.csv(results90245_df, "output/Seed_Zone2_Shifts_under_climate_change_2090_ssp245_KNN_method_extended.csv", row.names = FALSE)

#SSP90-585
results90585_df <- data.frame(
  Original_Longitude = df_pca$lon,
  Original_Latitude = df_pca$lat,
  Min_Distance = min_distances90585,
  Future_Longitude = Trns_grid_90585$lon[min_indices90585],
  Future_Latitude = Trns_grid_90585$lat[min_indices90585],
  Seed_Zone = mixcluster$color
)
write.csv(results90585_df, "output/Seed_Zone2_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv", row.names = FALSE)

############################################################################
#ssp245_70
two_zones_df <- read.csv("output/Seed_Zone2_Shifts_under_climate_change_2070_ssp245_KNN_method_extended.csv")
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

write.csv(arrow_data_70245,"output/arrow_data_Zone2_70245.csv")

# 获取中国省级边界数据
# China_provinces <- ne_states(country = "China", returnclass = "sf")
# taiwan_province <- ne_states(country = "taiwan", returnclass = "sf")
China_provinces <- st_read("中国_省.shp") #  本地导入  ne_states没有台湾省
arrow_data <- arrow_data_70245 

ggplot() +
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.1) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 自定义颜色映射
  scale_fill_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "black", linewidth = 0.4) +
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )

#ssp585_70
two_zones_df <- read.csv("output/Seed_Zone2_Shifts_under_climate_change_2070_ssp585_KNN_method_extended.csv")
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

write.csv(arrow_data_70585,"output/arrow_data_Zone2_70585.csv")

China_provinces <- st_read("中国_省.shp") #  本地导入shp  ne_states没有台湾省
arrow_data <- arrow_data_70585 

ggplot() +
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.1) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 自定义颜色映射
  scale_fill_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "black", linewidth = 0.4) +
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )

#ssp245_90
two_zones_df <- read.csv("output/Seed_Zone2_Shifts_under_climate_change_2090_ssp245_KNN_method_extended.csv")
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

write.csv(arrow_data_90245,"output/arrow_data_Zone2_90245.csv")

# 获取中国省级边界数据
# China_provinces <- ne_states(country = "China", returnclass = "sf")
# taiwan_province <- ne_states(country = "taiwan", returnclass = "sf")
China_provinces <- st_read("中国_省.shp") #  本地导入  ne_states没有台湾省
arrow_data <- arrow_data_90245 

ggplot() +
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.1) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 自定义颜色映射
  scale_fill_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "black", linewidth = 0.4) +
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )

#ssp585_90
two_zones_df <- read.csv("output/Seed_Zone2_Shifts_under_climate_change_2090_ssp585_KNN_method_extended.csv")
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

write.csv(arrow_data_90585,"output/arrow_data_Zone2_90585.csv")

China_provinces <- st_read("中国_省.shp") #  本地导入shp  ne_states没有台湾省
arrow_data <- arrow_data_90585 

ggplot() +
  # 绘制省级边界
  geom_sf(data = China_provinces, fill = "grey90", color = "black", lwd = 0.1) + 
  # 绘制种群区域的未来点
  geom_point(aes(x = two_zones_df$Future_Longitude, 
                 y = two_zones_df$Future_Latitude, 
                 col = as.factor(two_zones_df$Seed_Zone)), 
             pch = ".", cex = 0.1) +
  # 自定义颜色映射
  scale_color_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 自定义颜色映射
  scale_fill_manual(values = c(
    "1" = "#925372",
    "2" = "#AB8650"
  )) +
  # 绘制原始种群中心点
  geom_point(aes(x = centroids$Orig_Centroid_Lon, 
                 y = centroids$Orig_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在原始种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Orig_Centroid_Lon, 
                y = centroids$Orig_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "white") +
  # 绘制未来种群中心点
  geom_point(aes(x = centroids$Future_Centroid_Lon, 
                 y = centroids$Future_Centroid_Lat, 
                 fill = as.factor(centroids$Seed_Zone)), 
             pch = 21, size = 6, color = "black", stroke = 0.3) +
  # 在未来种群中心点的圆圈中添加数字
  geom_text(aes(x = centroids$Future_Centroid_Lon, 
                y = centroids$Future_Centroid_Lat, 
                label = centroids$Seed_Zone), 
            size = 3, color = "black") +
  # 绘制箭头连接原始和未来种群中心
  geom_segment(data = arrow_data, 
               aes(x = x, y = y, xend = xend, yend = yend), 
               arrow = arrow(length = unit(0.1, "cm")), 
               color = "black", linewidth = 0.4) +
  coord_sf(xlim = c(96, 123), ylim = c(18, 33), expand = FALSE) +
  xlab("Longitude") + 
  ylab("Latitude") +
  theme_bw() +
  theme(
    plot.background = element_blank(), 
    axis.text = element_text(size = 12,color = "black"), 
    axis.title = element_text(size = 14), 
    panel.grid.major = element_blank(), 
    legend.position = "none"
  )
