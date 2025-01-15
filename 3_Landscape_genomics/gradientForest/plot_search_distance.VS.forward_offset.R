# 加载必要的库
library(ggplot2)
library(dplyr)
library(stringr)

setwd("C:/Rstudio/RStudio/Workspace/gradientForest_2024_1m/Forward_Genetic_Offset_1m")

timerange <- "2081-2100"

# 定义文件路径列表
file_paths <- paste0(
  c(
    "Forward_Genetic_Offset_1m_10km_ssp245_",
    "Forward_Genetic_Offset_1m_10km_ssp585_",
    "Forward_Genetic_Offset_1m_20km_ssp245_",
    "Forward_Genetic_Offset_1m_20km_ssp585_",
    "Forward_Genetic_Offset_1m_50km_ssp245_",
    "Forward_Genetic_Offset_1m_50km_ssp585_",
    "Forward_Genetic_Offset_1m_100km_ssp245_",
    "Forward_Genetic_Offset_1m_100km_ssp585_",
    "Forward_Genetic_Offset_1m_200km_ssp245_",
    "Forward_Genetic_Offset_1m_200km_ssp585_",
    "Forward_Genetic_Offset_1m_500km_ssp245_",
    "Forward_Genetic_Offset_1m_500km_ssp585_"
  ),
  timerange, ".csv"
)
# 初始化空数据框
all_data <- data.frame()

# 遍历所有文件
for (file in file_paths) {
  # 使用 fread 读取文件
  df <- fread(file)
  
  # 提取文件名中的信息
  file_info <- strsplit(basename(file), "_")[[1]]
  search_distance <- gsub("km", "", file_info[4])  # 提取搜索距离
  search_distance <- ifelse(search_distance == "Unlimited", NA, as.numeric(search_distance))  # 处理 "Unlimited"
  
  scenario <- file_info[5]  # 提取情景
  time_range <- gsub(".csv", "", file_info[6])  # 提取时间范围
  
  # 添加到数据框
  df$search_distance <- search_distance
  df$scenario <- scenario
  df$time_range <- time_range
  
  # 合并数据
  all_data <- rbindlist(list(all_data, df), use.names = TRUE, fill = TRUE)
}

# 从 scenario 列提取搜索距离
all_data <- all_data %>%
  mutate(
    search_distance = str_extract(scenario, "\\d+km"),  # 提取包含数字和 "km" 的部分
    search_distance = as.numeric(gsub("km", "", search_distance))  # 去掉 "km" 并转换为数值
  )
# 计算每组的中位值和百分位数
summary_data <- all_data %>%
  group_by(search_distance, scenario, time_range) %>%
  summarise(
    median_offset = median(forwardOffset, na.rm = TRUE),
    lower_offset = quantile(forwardOffset, 0.25, na.rm = TRUE),
    upper_offset = quantile(forwardOffset, 0.75, na.rm = TRUE)
  )

ggplot(summary_data, aes(x = search_distance, y = median_offset, 
                         group = time_range, color = time_range)) +
  geom_line(linewidth = 0.75) +
  geom_ribbon(aes(ymin = lower_offset, ymax = upper_offset, fill = time_range), alpha = 0.2, color = NA) +
  geom_point(size = 2, shape = 21, fill = "white") +  # shape = 21 表示空心圆，fill 为内部填充色
  labs(
    title = NULL,
    x = "Search Distance (km)",
    y = "Forward Offset",
    color = NULL,
    fill = NULL
  ) +
  # 手动设置颜色
  scale_color_manual(values = c("#4B5CC4", "#9D2933")) +
  scale_fill_manual(values = c("#4B5CC4", "#9D2933")) +
  # 设置主题
  theme_bw() +
  theme(
    legend.position = c(0.85, 0.85),  # 图例放置在坐标区的右上角
    legend.background = element_rect(fill = "white", color = "gray50"),  # 图例背景白色，边框黑色
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 14, face = "bold"),
    # 移除内部网格线
    panel.grid.major = element_blank(),  # 移除主网格线
    panel.grid.minor = element_blank()   # 移除次网格线
  )
