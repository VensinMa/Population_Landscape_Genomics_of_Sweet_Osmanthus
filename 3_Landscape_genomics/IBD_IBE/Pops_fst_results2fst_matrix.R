# 读取 CSV 文件并生成 FST 矩阵
generate_fst_matrix <- function(input_file, output_file) {
  # 读取数据
  data <- read.csv(input_file, header=TRUE)

  # 提取唯一的组名
  unique_groups <- unique(c(data$Group1, data$Group2))

  # 创建一个全零矩阵，行列名为唯一的组名
  matrix <- matrix(0, nrow=length(unique_groups), ncol=length(unique_groups), dimnames=list(unique_groups, unique_groups))

  # 填充矩阵
  for (i in 1:nrow(data)) {
    row_name <- as.character(data[i, "Group1"])
    col_name <- as.character(data[i, "Group2"])
    matrix[row_name, col_name] <- as.numeric(data[i, "Mean_FST"])
    matrix[col_name, row_name] <- as.numeric(data[i, "Mean_FST"])  # 对称位置也填充
  }

  # 保存矩阵到文件
  write.csv(matrix, output_file)
}

# 使用示例
generate_fst_matrix("Pops_fst_results.csv", "Pops_fst_matrix.csv")
