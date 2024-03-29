library(raster)

# 读取包含经纬度信息的CSV文件
coordinates <- read.csv("E:/迅雷下载/MIROC6_30s/186sample_id_lat_lon_ele.csv")

# 设置TIFF文件目录
tif_directory <- "E:/迅雷下载/MIROC6_30s/future_climate"

# 创建新文件夹来存储结果文件
output_directory <- "extracted_data"
if (!dir.exists(output_directory)) {
  dir.create(output_directory, showWarnings = FALSE)
}

# 定义情景和时间范围的组合
scenario_time_combinations <- c(
  "ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
  "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100"
)

# 定义结果数据框的列顺序
column_order <- c(
  "lon", 
  "lat",
  paste("BIO", 1:19, sep = ""),
  "Elevation",
  "SRAD",
  "SOC",
  "PHH2O"
)

# 固定环境变量文件名称
fixed_env_files <- c(
  "srad_avg.tif" = "SRAD",
  "wc2.1_30s_elev.tif" = "Elevation",
  "phh2o_5cm_mean.tif" = "PHH2O",
  "soc_5cm_mean.tif" = "SOC"
)

# 循环遍历每个情景和时间范围的组合
for (combo in scenario_time_combinations) {
  # 创建一个数据框，包含经纬度列
  result <- data.frame(lon = coordinates$lon, lat = coordinates$lat)
  
  # 先处理固定环境变量
  for (file_name in names(fixed_env_files)) {
    file_path <- file.path(tif_directory, file_name)
    if (file.exists(file_path)) {
      env_data <- raster(file_path)
      values <- extract(env_data, coordinates[, c("lon", "lat")])
      result[[fixed_env_files[file_name]]] <- values
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  # 处理生物气候变量
  for (i in 1:19) {
    env_name <- paste("BIO", i, sep = "")
    file_name <- sprintf("wc2.1_30s_bioc_MIROC6_%s_%s.tif", combo, env_name)
    file_path <- file.path(tif_directory, file_name)
    
    if (file.exists(file_path)) {
      env_data <- raster(file_path)
      values <- extract(env_data, coordinates[, c("lon", "lat")])
      result[[env_name]] <- values
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  # 调整结果数据框的列顺序确保一致性
  result <- result[, column_order]
  
  # 构造输出文件路径
  output_file_name <- paste("future_climate", combo, "O.fragrans.csv", sep = "_")
  output_file_path <- file.path(output_directory, output_file_name)
  
  # 导出结果数据框为 CSV 文件
  write.csv(result, output_file_path, row.names = FALSE)
}
