# 加载必要的库
library(ggplot2)
library(dplyr)
library(ggpubr)  # 提供 stat_compare_means 函数
library(ggforce)  # 提供 geom_sina 和 geom_half_violin
library(gghalves)

# 1. 加载主数据：这里的数据就像宇宙的起点，是一切分析的开端
data <- read.table("C:/Users/mawenxin/workspace/est-sfs/indv_GT_stats_res.197samples_filtered_3_outgroup.polarized.snpeff_LOF.txt", 
                   header = TRUE, stringsAsFactors = FALSE)  # 设置 header 以读取列名

# 2. 加载分组文件：没有分组，数据也就失去了灵魂
group_data <- read.table("C:/Users/mawenxin/workspace/est-sfs/prepare_est-sfs/200sample_4lineage.id", 
                         header = FALSE, stringsAsFactors = FALSE)
colnames(group_data) <- c("indv", "lineage")  # 为列添加标签，indv 是个体，lineage 是群体

# 3. 合并主数据和分组信息：数据整合，像星云汇聚成银河
merged_data <- merge(data, group_data, by = "indv")

# 4. 移除 Outgroup 分组
merged_data <- merged_data %>% filter(lineage != "Outgroup")  # 使用 filter 移除 Outgroup

# 5. 计算 hom1 / total 的比值：一个简单的数学操作，但能揭示深远的意义
merged_data$p11 <- merged_data$hom1 / 808574  # 分母是总值，固定数值确保比例一致

# 6. 按 lineage 计算均值并排序：让数据按照贡献大小排排坐，吃果果
lineage_means <- merged_data %>%
  group_by(lineage) %>%
  summarise(mean_p11 = mean(p11, na.rm = TRUE)) %>%  # 忽略 NA 值，计算平均值
  arrange(mean_p11)  # 按均值升序排列

# 7. 将 lineage 转为因子，按照均值排序：排序有序，数据就像军队一样整齐
merged_data$lineage <- factor(merged_data$lineage, levels = lineage_means$lineage)

# 按 lineage 手动指定顺序
merged_data$lineage <- factor(
  merged_data$lineage, 
  levels = c("West-YN", "West-GZ", "Centrol", "East")  # 指定分组顺序
)

# 检查因子顺序是否正确
levels(merged_data$lineage)

# 8. 定义颜色（根据分组名称设置）：色彩不只是美学的体现，也是数据的指路牌
fill_colors <- c(
  "East" = "#E5A84B",      # 黄河琉璃：温暖的亮黄色
  "Centrol" = "#DD6B4F",   # 佰赭：稳重的赭石红
  "West-GZ" = "#509296",   # 扁青：平静的青绿色
  "West-YN" = "#DD7694"    # 苏梅：带有雨水的淡粉红
)

# 9. 绘制箱线图，加入显著性差异分析：在数据世界里，用图表讲故事
p <- ggplot(merged_data, aes(x = lineage, y = p11, fill = lineage)) +
  geom_boxplot(colour = "black",alpha = 0.8) +  # 使用箱线图展示数据分布
  scale_fill_manual(values = fill_colors) +  # 应用自定义颜色
  labs(
    title = "p11 by Lineage (Excluding Outgroup)",  # 图表标题
    #x = "Lineage",  # x 轴标签
    y = "LOF"  # y 轴标签
  ) +
  theme_bw() +  # 设置简洁主题，让数据更加突出
  theme(colour = "black",
        axis.text = element_text(colour = "black"),  # 设置坐标轴文字颜色为黑色
        axis.title = element_text(colour = "black"),  # 设置坐标轴标题颜色为黑色
        legend.position = "none") + # 隐藏图例，减少干扰  # x 轴标签旋转，防止拥挤
  #theme( axis.text.x = element_text(angle = 45, hjust = 1)) +  # x 轴标签旋转，防止拥挤
  stat_compare_means(
    comparisons = list(
      c("Centrol", "East"), 
      c("West-GZ", "Centrol"),
      c("West-GZ", "East"),
      c("West-YN", "West-GZ"),
      c("West-YN", "Centrol"),
      c("West-YN", "East")
    ),
    label = "p.signif",  # 显示显著性差异的星号标记
    method = "t.test"  # 使用 t 检验
  )


# 9. 绘制一半小提琴图 + 箱线图 + 散点图，加入显著性差异分析
p <- ggplot(merged_data, aes(x = lineage, y = p11, fill = lineage)) +
  # 添加右半小提琴图
  geom_half_violin(
    side = "r",  # 右侧小提琴图
    position = position_nudge(x = .15),  # 向右微调
    alpha = 0.9,  # 设置透明度
    colour = NA,  # 去除边框线
    scale = "area",  # 小提琴宽度按数据分布调整
    width = 1  # 加宽小提琴图
  ) +
  # 添加箱线图
  geom_boxplot(
    width = 0.07,
    position = position_nudge(x = .15),  # 对齐小提琴图
    fill = "white",  # 箱线图内部为白色
    colour = "black",  # 箱线边框为黑色
    outlier.shape = NA,  # 去除异常点显示
    size = 0.21  # 调整箱线图边框线条粗细
  ) +
  # 添加散点图，结合抖动和位置偏移
  geom_jitter(
    aes(color = lineage), 
    size = 1, 
    position = position_jitterdodge(jitter.width = 0.4, dodge.width = 0.5),  # 调整抖动和偏移宽度
    alpha = 0.8
  ) +
  # 添加显著性差异标记
  stat_compare_means(
    comparisons = list(
      c("Centrol", "East"), 
      c("West-GZ", "Centrol"),
      c("West-GZ", "East"),
      c("West-YN", "West-GZ"),
      c("West-YN", "Centrol"),
      c("West-YN", "East")
    ),
    label = "p.signif",  # 显示星号表示显著性
    method = "t.test",  # 使用 t 检验
    step.increase = 0.05  # 调整显著性标记的高度
  ) +
  # 定义填充颜色
  scale_fill_manual(values = fill_colors) +
  # 定义散点图颜色
  scale_color_manual(values = fill_colors, guide = "none") +  # 隐藏散点图图例
  # 设置图形主题
  theme_bw(base_size = 10) +
  labs(
    title = expression("Homozygous genotypes (" * italic(p)[1/1] * ")"), # 图标题
    x = NULL,  # 去掉 x 轴标题
    y = "LOF"  # y 轴标题
  ) +
  #scale_y_continuous(limits = c(0, max(merged_data$p11))) +  # 调整 y 轴范围，可以不显示显著性部分
  theme(
    axis.title = element_text(size = 10, face = "bold"),  # 坐标轴标题样式
    axis.text = element_text(size = 10, color = "black", face = "bold"),  # 坐标轴刻度文字样式
    #axis.text.x = element_text(angle = 45, hjust = 1, face = "bold.italic"),  # x 轴刻度倾斜样式
    legend.position = "none"  # 隐藏图例
  )


# 10. 显示图形：让数据以优雅的姿态亮相
print(p)

# 11. 保存图形为图片（可选）：数据成果需要存档，以备分享和膜拜
pdf(file = "LOF_p11_4LIN_3_outgroup.pdf", height = 6, width = 4,)
print(p)
dev.off()