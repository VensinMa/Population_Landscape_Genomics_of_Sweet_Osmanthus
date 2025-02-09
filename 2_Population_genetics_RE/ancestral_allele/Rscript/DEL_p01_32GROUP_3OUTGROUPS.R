# 加载必要的库
library(dplyr)  # 用于数据处理
library(ggplot2)  # 用于数据可视化
library(ggpubr)  # 提供显著性差异分析功能
library(RColorBrewer)  # 用于生成调色板

# 1. 加载主数据：读取包含样本和统计信息的数据文件
data <- read.table("C:/Users/mawenxin/workspace/sift-lyjg/indv_GT_stats_res.deleterious_197.txt", 
                   header = TRUE, stringsAsFactors = FALSE)  # 设置 header 为 TRUE 以读取列名

# 2. 提取群体信息：从 indv 列中提取 "-" 前的字母作为群体标识
data$group <- sub("-.*", "", data$indv)  # 使用正则表达式截取字符串，"-.*" 表示匹配 "-" 后的所有内容

# 3. 过滤掉外类群 “O”
data <- data %>% filter(group != "O")  # 移除 group 列中等于 "O" 的行

# 4. 计算比例：将 het 除以 total，得到每个样本的比例
data$p01 <- data$het / 808574  # 直接按列进行数学运算

# 5. 按群体计算均值：分组计算每个群体的 p01 均值
group_means <- data %>%
  group_by(group) %>%  # 按 group 分组
  summarise(mean_p01 = mean(p01, na.rm = TRUE))  # 计算每组的平均值，na.rm = TRUE 忽略 NA 值

# 6. 按均值排序群体：根据 mean_p01 排序群体的顺序
group_means <- group_means %>%
  arrange(mean_p01)  # 按 mean_p01 升序排列

# 7. 因子化群体：将群体列转为因子并按照均值排序
data$group <- factor(data$group, levels = group_means$group)  # 设置因子的 levels，按照均值排序后的顺序

# 8. 绘制箱线图
data$lineage <- "East"  # 默认分组为 East
data$lineage[data$group %in% c("DRS", "LCJ")] <- "West-YN"
data$lineage[data$group %in% c("RX", "JMX", "ZJS", "HYX", "HJLL")] <- "West-GZ"
data$lineage[data$group %in% c("DA", "DST", "EJ", "ZJP")] <- "Centrol"

# 自定义填充颜色
custom_colors <- c(
  "East" = "#E5A84B",      # 黄河琉璃：温暖的亮黄色
  "Centrol" = "#DD6B4F",   # 佰赭：稳重的赭石红
  "West-GZ" = "#509296",   # 扁青：平静的青绿色
  "West-YN" = "#DD7694"    # 苏梅：带有雨水的淡粉红
)

# 绘制 箱线图 + 散点图
p = ggplot(data, aes(x = group, y = p01, fill = lineage)) +
  # 添加箱线图
  geom_boxplot(
    width = 0.8,  # 箱线图宽度
    colour = "black",  # 箱线边框为黑色
    outlier.shape = NA,  # 去除异常点显示
    size = 0.15  # 调整箱线图边框线条粗细
  ) +
  # 添加散点图，结合抖动和位置偏移
  geom_jitter(
    aes(fill = lineage),  # 填充颜色映射到 lineage
    shape = 21,  # 点形状，允许独立设置边框和填充
    color = "black",  # 边框颜色设为黑色
    size = 1,  # 点大小
    stroke = 0.1,  # 边框线条粗细
    position = position_jitterdodge(
      jitter.width = 0.2, 
      dodge.width = 0.15
    ),  
    alpha = 0.8  # 设置透明度
  ) +
  # 定义填充颜色
  scale_fill_manual(values = custom_colors) +
  # 定义散点图颜色
  scale_color_manual(values = custom_colors, guide = "none") +  # 隐藏散点图图例
  # 设置图形主题
  theme_bw(base_size = 10) +
  labs(
    title = expression("Heterozygous genotypes (" * italic(p)[0/1] * ")"), # 图标题
    x = NULL,  # x 轴标题
    y = "DEL"  # y 轴标题
  ) +
  theme(
    axis.title = element_text(size = 10, face = "bold"),  # 坐标轴标题样式
    axis.text = element_text(size = 10, color = "black", face = "bold"),  # 坐标轴刻度文字样式
    axis.text.x = element_text(angle = 45, hjust = 1),  # x 轴刻度倾斜
    legend.position = "none"  # 隐藏图例
  )

# 9. 显示图形
print(p)

# 10. 保存图片
pdf(file = "DEL_p01_32GROUP_3OUTGROUPS.pdf", height = 5, width = 7.5,)
print(p)
dev.off()