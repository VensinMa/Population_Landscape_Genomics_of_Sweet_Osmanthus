# 加载raster库
library(raster)

# 指定包含未来气候数据的文件夹路径
folder_path <- "E:/迅雷下载/MIROC6_30s/"

# 获取文件夹中所有后缀为.tif的文件名
tif_files <- list.files(folder_path, pattern = "\\.tif$", full.names = TRUE)

# 遍历每个.tif文件
for (file in tif_files) {
  # 读取文件为raster brick对象
  dataset <- brick(file)
  
  # 为当前文件创建一个新的文件夹来存储BIO图层
  scene_folder <- file.path(folder_path, tools::file_path_sans_ext(basename(file)))
  if(!dir.exists(scene_folder)) {
    dir.create(scene_folder)
  }
  
  # 确保文件包含预期的19个BIO图层
  band_count <- nlayers(dataset)
  if(band_count == 19) {
    # 提取每个BIO图层
    for (i in 1:band_count) {
      # 读取当前BIO图层数据
      band_data <- subset(dataset, i)
      
      # 构建输出文件名，包含场景信息和BIO编号
      output_file_name <- paste0(tools::file_path_sans_ext(basename(file)), "_BIO", i, ".tif")
      output_file_path <- file.path(scene_folder, output_file_name)
      
      # 写入数据，将当前BIO图层保存为独立的.tif文件
      writeRaster(band_data, filename = output_file_path, format = "GTiff", overwrite = TRUE)
    }
  } else {
    cat("Warning: File", basename(file), "does not contain 19 layers. It has", band_count, "layers.\n")
  }
}
