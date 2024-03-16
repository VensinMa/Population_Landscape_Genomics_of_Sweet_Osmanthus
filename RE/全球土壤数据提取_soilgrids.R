###########################    加载必要的库  #################################
if (!requireNamespace("geodata", quietly = TRUE)) {
  install.packages("geodata")
}

if (!requireNamespace("raster", quietly = TRUE)) {
  install.packages("raster")
}

if (!requireNamespace("sp", quietly = TRUE)) {
  install.packages("sp")
}

if (!requireNamespace("rgdal", quietly = TRUE)) {
  install.packages("rgdal")
}

library(geodata) # 用于获取全球土壤数据
library(raster)  # 用于处理栅格数据
library(sp)      # 用于空间数据的类和方法
library(rgdal)   # 用于地理数据的读写和转换

##############################  设置工作目录  ################################
# 设置工作目录的文件夹路径
folder_path <- "soilgrids"

# 检查文件夹是否存在
if (!dir.exists(folder_path)) {
  # 如果文件夹不存在，则创建文件夹
  dir.create(folder_path, recursive = TRUE)
}

# 设置工作目录
setwd(folder_path)

## https://www.soilgrids.org/
################################################################################

# 获取全球土壤pH值数据
phh2o <- soil_world(var="phh2o", depth=5, stat="mean", path=tempdir())
# 可视化土壤pH值数据
plot(phh2o, ylim = c(-60,90), xlim = c(-180,180),
     col = colorRampPalette(c("#3288BD","#66C2A5","#FEE08B","#D53E4F"))(250))

# 将土壤pH值数据转换为RasterLayer对象
phphh2o_raster <- raster(phh2o)
plot(phphh2o_raster)

# 导出RasterLayer为TIF格式 
writeRaster(phphh2o_raster, filename="phh2o_5cm_mean.tif", format="GTiff", overwrite=TRUE)

# 读取地理位置数据
data = read.csv("186sample_id_lat_lon_ele.csv")
coordinates(data) = c("lon", "lat") # 设置数据框的坐标

# 提取地理位置对应的土壤pH值
sgphh2o <- extract(x = phphh2o_raster, y = data)
datasgphh2o <- data.frame (data,sgphh2o) # 合并地理位置数据与土壤pH值
write.csv(datasgphh2o,"phh2o_5cm_mean.csv") # 保存pH值数据到CSV文件

###############################################################################

# 获取全球土壤有机碳（SOC）数据
soc <- soil_world(var="soc", depth=5, stat = "mean", path=tempdir())
# 可视化SOC数据
plot(soc, ylim = c(-60,90), xlim = c(-180,180),
     col = colorRampPalette(c("#3288BD","#66C2A5","#FEE08B","#D53E4F"))(250))

# 将SOC数据转换为RasterLayer对象
soc_raster <- raster(soc)
plot(soc_raster)

# 导出RasterLayer为TIF格式 
writeRaster(soc_raster, filename="soc_5cm_mean.tif", format="GTiff", overwrite=TRUE)

# 提取地理位置对应的SOC数据
sgsoc <- extract(x = soc_raster, y = data)
datasgsoc <- data.frame (data, sgsoc) # 合并地理位置数据与SOC数据
# 保存SOC数据到另一个CSV文件，避免覆盖pH值数据
write.csv(datasgsoc,"soc_5cm_mean.csv") 

