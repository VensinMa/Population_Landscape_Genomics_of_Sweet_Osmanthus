library(terra)

# 设置工作目录
getwd()
setwd("C:/Rstudio/RStudio/Workspace/GEA_2024/lfmm")
getwd()

# 读取CSV文件，假设第一列为样品名称，后续为经纬度等信息
coordinates <- read.csv("E:/OneDrive/文档/194sample_lon_lat_ele.csv")
str(coordinates)
## 'data.frame':	194 obs. of  4 variables:
##  $ ind: chr  "CP-1" "CP-2" "CP-3" "CP-4" ...
##  $ lon: num  118 118 118 118 118 ...
##  $ lat: num  28 28 28 28 28 ...
##  $ ele: num  418 442 305 273 270 270 262 255 252 238 ...

# 不同时期及场景下 BIO1-BIO19 环境图层
tif_directory <- "D:/DATA/worldclim/Future climate data/MIROC6_30s/future_climate"

# 提取样品名称列（假设第一列为样品名称）
sample_ids <- coordinates[, 1]

# 创建输出目录
output_directory <- "extracted_data"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, showWarnings = FALSE)
}

# 定义情景和时间范围
scenario_time_combinations <- c("current",
     "ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
     "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100"
)

# 定义列顺序
column_order <- c(
  "ind",  # 将样品名称作为第一列
  "lon", 
  "lat",
  paste("BIO", 1:19, sep = ""),
  "Elevation",
  "SRAD",
  "SOC",
  "PHH2O"
)

# 固定环境变量文件
fixed_env_files <- list(
  SRAD = "srad_avg.tif",
  Elevation = "wc2.1_30s_elev.tif",
  PHH2O = "phh2o_5cm_mean.tif",
  SOC = "soc_5cm_mean.tif"
)

# 提取栅格数据的函数
extract_raster_values <- function(file_path, coords) {
  if (file.exists(file_path)) {
    rast_data <- rast(file_path)
    return(extract(rast_data, coords[, c("lon", "lat")])[,2])
  } else {
    warning(paste("File not found:", file_path))
    return(rep(NA, nrow(coords)))
  }
}

# 循环处理每个情景
for (combo in scenario_time_combinations) {
  # 创建结果数据框，首先添加样品名称和经纬度信息
  result <- data.frame(ind = sample_ids, lon = coordinates$lon, lat = coordinates$lat)
  
  # 提取固定环境变量
  for (var_name in names(fixed_env_files)) {
    file_path <- file.path(tif_directory, fixed_env_files[[var_name]])
    result[[var_name]] <- extract_raster_values(file_path, coordinates)
  }
  
  # 提取生物气候变量
  for (i in 1:19) {
    bio_name <- paste("BIO", i, sep = "")
    file_name <- sprintf("wc2.1_30s_bioc_MIROC6_%s_%s.tif", combo, bio_name)
    file_path <- file.path(tif_directory, file_name)
    result[[bio_name]] <- extract_raster_values(file_path, coordinates)
  }
  
  # 调整列顺序
  result <- result[, column_order]
  
  # 保存结果
  output_file_name <- sprintf("Climate_%s_194samples.csv", combo)
  write.csv(result, file.path(output_directory, output_file_name), row.names = FALSE)
}
