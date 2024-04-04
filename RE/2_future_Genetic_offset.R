####################     
# 定义时期列表
periods <- c("ssp245_2041-2060", "ssp245_2061-2080", "ssp245_2081-2100",
             "ssp585_2041-2060", "ssp585_2061-2080", "ssp585_2081-2100")

# 确保存放结果的目录存在
if (!dir.exists("Genetic_Offset")) {
  dir.create("Genetic_Offset")
}

# 循环遍历每个时期
for (period in periods) {
  # 读取未来气候数据
  file_name <- paste0("extracted_future_data/future_climate_", period, "_O.fragrans.csv")
  future_data <- read.csv(file_name)
  
  # 筛选数据，并删除包含NA的行
  future_data <- future_data[complete.cases(future_data[, 3:25]), 1:25]
  
  # 使用梯度森林模型进行环境梯度预测
  future_data_pred <- cbind(future_data[,c("lon","lat")], 
                            predict(gf.mod, future_data[,3:25]))
  
  # 计算遗传偏移
  genOffsetAll <- sqrt(rowSums((future_data_pred[, 3:ncol(future_data_pred)] - 
                                  all_tgrid[, 3:ncol(all_tgrid)])^2))
  
  # 合并遗传偏移值到坐标数据中
  Offset <- cbind(future_data_pred[,c("lon","lat")], genOffsetAll)
  colnames(Offset)[3] <-"offset"
  
  # 保存遗传偏移结果为CSV文件
  output_file_name <- paste0("Genetic_Offset/", period, "_genetic_offset.csv")
  write.csv(Offset, output_file_name, quote=FALSE, row.names=FALSE)
}
