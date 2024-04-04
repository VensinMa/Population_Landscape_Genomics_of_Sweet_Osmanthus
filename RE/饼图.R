# 加载ggplot2包
library(ggplot2)
library(dplyr)
library(ggsci)
# 数据集1
data_frame1 <- data.frame(
  category = c("downstream", "intronic", "intergenic",  "splicing", "upstream", "exonic","upstream;downstream"),
  value = c(1102763, 592470, 7548615, 1906150, 3022, 1247057, 151190)
)

# 计算百分比
data_frame1 <- data_frame1 %>%
  mutate(percent = value / sum(value) * 100,
         label = paste0(category, "\n", value, "\n", round(percent, 2), "%"))

library(cols4all)
c4a_gui()

colors = c4a("paired", 9)
colors = c("#1F78B4", "#FB9A99", "#33A02C", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6")
# 绘制数据集1的饼图并添加标签
ggplot(data_frame1, aes(x = "", y = value, fill = category)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.title = element_blank()) +
  geom_label(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(fill = "Category", title = "Variants by Genomic Region Pie Chart")+
  scale_fill_manual(values = colors)


# 数据集2
data_frame2 <- data.frame(
  category = c("nonsynonymous", "stopgain", "stoploss", "synonymous", "unknown"),
  value = c(322622, 6955, 769, 261796, 328)
)

# 计算百分比
data_frame2 <- data_frame2 %>%
  mutate(percent = value / sum(value) * 100,
         label = paste0(category, "\n", value, "\n", round(percent, 2), "%"))

pastel1
colors2 = c4a("brewer.pastel1", 8)
# 绘制数据集2的饼图并添加标签
ggplot(data_frame2, aes(x = "", y = value, fill = category)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar("y", start = 0) +
  theme_void() +
  theme(legend.title = element_blank()) +
  geom_label(aes(label = label), position = position_stack(vjust = 0.5)) +
  labs(fill = "Category", title = "Exonic Variant Functions Pie Chart")+
  scale_fill_manual(values = colors2)
