# 加载必要的库
library(ggplot2)

# 定义五种基准颜色
colors <- c("#00ff48", "#0000FF", "#ff0000", "#FFFF00", "#00ff48")
colors <- c("#0EB83A", "#0034DE", "#FF461F", "#FAFF14", "#0EB83A")
colors <- c("#698E6A", "#69549D", "#B15A43", "#FAFF72", "#698E6A")

# 使用 colorRampPalette 生成60个渐变色
gradient_colors <- colorRampPalette(colors)(60)
gradient_colors <- colorRampPalette(colors)(120)
# 打印生成的渐变色
print(gradient_colors)

# 可视化渐变色
barplot(rep(1, 60), col = gradient_colors, border = NA, space = 0,
        main = "60-step Gradient Color Bar", xlab = "Color Index")

barplot(rep(1, 120), col = gradient_colors, border = NA, space = 0,
        main = "120-step Gradient Color Bar", xlab = "Color Index")
# 加载生成的 bearing_6_degree_counts.csv 数据
data <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024_1m/Forward_Genetic_Offset_1m/bearing_6_degree_counts.csv")
data <- read.csv("C:/Rstudio/RStudio/Workspace/gradientForest_2024_1m/Forward_Genetic_Offset_1m/bearing_3_degree_counts.csv")


# 计算方位角区间的中值
bin_interval <- 3  # 设置区间间隔
# 计算 bearing_group_mid 直接通过 bearing_bin 起点加上半个区间
data$bearing_group_mid <- data$bearing_bin + (bin_interval / 2)

# 将生成的渐变色映射到 bearing_group_mid 的分组
data$color <- gradient_colors[cut(data$bearing_group_mid, breaks = seq(0, 360, length.out = 61), labels = FALSE)]
data$color <- gradient_colors[cut(data$bearing_group_mid, breaks = seq(0, 360, length.out = 121), labels = FALSE)]

# 绘制极坐标直方图
ggplot(data, aes(x = bearing_group_mid, y = count, fill = color)) +
  geom_bar(stat = "identity", width = bin_interval) +  # 根据区间间隔调整条宽
  coord_polar(start = 0) +
  scale_x_continuous(breaks = c(0, 90, 180, 270), 
                     labels = c("0°", "90°", "180°", "270°")) +
  labs(title = "Bearing Group Counts (Polar Histogram)",
       x = "Bearing Group (Degrees)",
       y = "Counts") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 12, face = "bold"),
        legend.position = "none") +  # 隐藏图例
  scale_fill_identity()  # 保持颜色不变
