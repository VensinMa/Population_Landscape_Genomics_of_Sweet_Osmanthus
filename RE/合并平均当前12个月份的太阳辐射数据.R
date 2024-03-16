library(raster)
library(dplyr)

# 定义文件路径
file_paths <- list.files(path="D:/DATA/worldclim/Historical climate data/wc2.1_30s_srad", 
                         pattern="srad_\\d+.tif$", full.names=TRUE)

# 读取所有TIFF文件为RasterLayer对象，并存储在一个列表中
srad_layers <- lapply(file_paths, raster)

# 计算这些RasterLayer对象的平均值
srad_avg <- stack(srad_layers) %>% calc(fun=mean)

# 可视化平均值
plot(srad_avg)

# 保存平均值栅格数据为新的TIFF文件
writeRaster(srad_avg, filename="D:/DATA/worldclim/Historical climate data/wc2.1_30s_srad/srad_avg.tif", format="GTiff", overwrite=TRUE)
