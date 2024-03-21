########################  1 从多波段tif提取各波段气候数据 ######################
# 加载raster库，用于处理地理信息数据
library(raster)

# 指定文件夹路径，为包含未来气候数据的文件夹路径
folder_path <- "D:/DATA/worldclim/Future climate data/CMIP6_0.5"

# pattern = "//.tif$" 用于匹配文件夹中所有后缀为.tif的文件
tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# 遍历每个.tif文件
for (file in tif_files) {
  # 创建以.tif文件名命名的文件夹，用于存储各个波段的数据
  output_folder <- file.path(folder_path, tools::file_path_sans_ext(basename(file)))
  dir.create(output_folder, showWarnings = FALSE)
  # 打开.tif文件作为brick对象，brick是用于处理多波段栅格数据的对象
  dataset <- brick(file)
  # 获取波段数量
  band_count <- nlayers(dataset)
  # 提取每个波段的数据并保存在相应的文件夹中
  for (i in 1:band_count) {
    # 读取波段数据，subset函数用于提取指定波段
    band_data <- subset(dataset, i)
    # 构建输出文件名，每个波段将保存为单独的.tif文件
    output_file <- file.path(output_folder, paste0("bio", i, ".tif"))
    # 写入数据，将每个波段的数据保存为独立的.tif文件
    writeRaster(band_data, filename = output_file, format = "GTiff", overwrite = TRUE)
  }
}

################ 2 从19个气候因子图层中提取指定经纬度数据 ######################
## 未来的多波段图层经提取后，可和当前bio1-bio19一样提取
# 加载raster库，用于处理地理信息数据
library(raster)
getwd()
#  "C:/RStudio/RStudio/Workspace/Biomod2"
# 读取包含经纬度信息的CSV文件
coordinates <- read.csv("O.fragrans_cunzai.csv") 
dim(coordinates)
# 设置TIFF文件目录
tif_directory <- "D:/DATA/worldclim/Historical climate data/wc2.1_30s_bio"
# 指定要提取的TIFF文件名列表
specified_files <- c(paste("wc2.1_30s_bio_", 1:19, ".tif", sep = ""),
                     "wc2.1_30s_elev.tif", 
                     "wc2.1_30s_Aspect.tif",
                     "wc2.1_30s_Slope.tif")
# 映射的列名
specified_cols  =  c(paste("bio_", 1:19, sep = ""),
                     "Elev", "Aspect", "Slope" )
# 创建一个命名向量，将文件名映射到列名
column_mapping <- setNames(
  specified_cols,
  specified_files
)

# 创建一个数据框，包含经纬度列
result <- data.frame(
  lon = coordinates$lon, 
  lat = coordinates$lat
)

# 循环遍历指定的 TIFF 文件
for (tif_file in specified_files) {
  # 从映射表中获取列名
  env_name <- column_mapping[tif_file]
  
  # 读取环境数据 TIFF 文件
  env_data <- raster(file.path(tif_directory, tif_file))
  
  # 提取经纬度位置的值并命名列为环境数据的名称
  values <- extract(env_data, coordinates[, c("lon", "lat")])
  
  # 将提取的值添加到结果数据框，使用 env_name 作为列名
  result[[env_name]] <- values
}

# 导出结果数据框为 CSV 文件
write.csv(result, "wc2.1_30s_bioc_present_O.fragrans.csv", row.names = FALSE)
