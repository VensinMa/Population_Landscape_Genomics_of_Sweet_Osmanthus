library(ggplot2)
library(vegan)  # 确保你已经安装了vegan包

# 假设of.rda是你的RDA分析结果
rda_scores <- scores(of.rda, display = "sites", scaling = 1)  # 获取样本得分
species_scores <- scores(of.rda, display = "species", scaling = 1)  # 获取SNP得分
bp_scores <- scores(of.rda, display = "bp", scaling = 1)  # 获取预测变量得分

# 将RDA结果转换为数据框
sample_data <- data.frame(rda_scores, group = env_data$group)
bp_data <- data.frame(bp_scores, Label = rownames(bp_scores))

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

print(p)


library(ggplot2)
library(ggrepel)
library(vegan)  # 确保vegan包已安装

ggplot() +
  # 绘制样本点
  geom_point(data = sample_data, aes(x = RDA1, y = RDA2, fill = group), 
             shape = 21, size = 2, color = "gray20", stroke = 0.0001) +
  scale_fill_manual(values = color_mapping) +  # 设置颜色映射
  # 绘制预测变量向量
  geom_segment(data = bp_data, aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               arrow = arrow(length = unit(0.1, "cm")), color = "#0868ac") +
  # 使用 geom_text_repel 防止文本重叠
  geom_text_repel(data = bp_data, aes(x = RDA1, y = RDA2, label = Label), 
                  color = "black", size = 2, max.overlaps = 30) +
  labs(x = "RDA1", y = "RDA2", title = NULL) +
  theme_bw() +
  theme(legend.position = "bottom")
