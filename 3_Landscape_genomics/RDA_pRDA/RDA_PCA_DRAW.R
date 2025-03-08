library(ggplot2)
library(vegan)  # 确保你已经安装了vegan包
load("of.rda.1002.RData")
env_data = read.csv("Climate_current_311samples.csv",
                    header = T, row.names = 1)

# 假设of.rda是你的RDA分析结果
rda_scores <- scores(of.rda, display = "sites", scaling = 1)  # 获取样本得分
species_scores <- scores(of.rda, display = "species", scaling = 1)  # 获取SNP得分
bp_scores <- scores(of.rda, display = "bp", scaling = 1)  # 获取预测变量得分

# 将RDA结果转换为数据框
sample_data <- data.frame(rda_scores, group = env_data$group)
bp_data <- data.frame(bp_scores, Label = rownames(bp_scores))

# 假定colors包含了三个颜色代码
colors <- c("#FDD379", "#E26472", "#A6DAEF")

# 设置env_data$pop的水平
levels(env_data$group) <- c("North", "Central", "South")

# 创建一个命名向量，将分组名称映射到颜色
color_mapping <- setNames(colors, levels(env_data$group))

# 将env_data$group中的每个分组名称替换为对应的颜色
env_data$group_color <- color_mapping[env_data$group]

# 打印结果，查看分组对应的颜色
bg = env_data$group_color


# ggplot作图
ggplot() +
  # 绘制样本点
  geom_point(data = sample_data, aes(x = RDA1, y = RDA2, fill = group), 
             shape = 21, size = 1.5, color = "gray20",stroke = 0.0001) +
  scale_fill_manual(values = color_mapping) +  # 设置颜色映射
  # 绘制预测变量向量
  geom_segment(data = bp_data, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.2, "cm")), color = "#0868ac") +
  geom_text(data = bp_data, aes(x = RDA1, y = RDA2, label = Label), 
            color = "#0868ac", size = 2, hjust = 0.5, vjust = -0.5) +
  labs(x = "RDA1", y = "RDA2", title = "RDA Analysis") +
  theme_bw() +
  theme(legend.position = "bottom")

