# 加载必要的库
library(sp)
library(raster)

# 读取地理位置数据
data <- read.csv("C:/RStudio/RStudio/Workspace/GEA/soilgrids/186sample_id_lat_lon_ele.csv")

# 将数据框转换为SpatialPointsDataFrame对象
coordinates(data) <- ~lon+lat

# 定义太阳辐射数据的文件路径模板
srad_path_template <- "D:/DATA/worldclim/Historical climate data/wc2.1_30s_srad/wc2.1_30s_srad_%02d.tif"

# 初始化一个向量来存储各个月份的太阳辐射值
srad_values <- vector("list", 12)

# 循环读取每个月的太阳辐射数据并提取指定地点的值
for (i in 1:12) {
  # 生成当前月份的文件路径
  srad_path <- sprintf(srad_path_template, i)
  # 读取太阳辐射数据
  srad_raster <- raster(srad_path)
  # 提取地点的太阳辐射值
  srad_values[[i]] <- extract(srad_raster, data)
}

# 将提取的太阳辐射值列表转换为矩阵
srad_matrix <- do.call(cbind, srad_values)

# 计算每一行（即每个地点）的平均太阳辐射值
srad_mean <- rowMeans(srad_matrix, na.rm = TRUE)

# 将平均值添加到原始数据框中
data$srad_mean <- srad_mean

# 查看更新后的数据框头部
head(data)

# 将更新后的数据框保存为新的CSV文件
write.csv(data, "C:/RStudio/RStudio/Workspace/GEA/soilgrids/186sample_id_lat_lon_ele_srad_mean.csv", row.names = FALSE)
