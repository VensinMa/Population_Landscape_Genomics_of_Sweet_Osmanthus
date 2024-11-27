## Rscript Pops_mean_env_vars.R

# 文件名
input_file <- "194sample_id_group_lat_lon_ele_with_env_vars.csv"  # 输入文件
output_file <- "32pop_means_env_vars.csv"  # 输出文件

# 加载数据
data <- read.csv(input_file, header = TRUE)

# 按照 pop 分组求平均值（忽略非数值列）
grouped_data <- aggregate(. ~ pop, data = data[, -1], FUN = mean)

# 保存分组平均值到文件
write.csv(grouped_data, output_file, row.names = FALSE)

# 提示完成
cat("按群体计算环境变量平均值，结果已保存到文件：", output_file, "\n")
