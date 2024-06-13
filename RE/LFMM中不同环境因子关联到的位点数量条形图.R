# 安装并加载ggplot2包

library(ggplot2)

# 数据
data <- data.frame(
  Variable = c("BIO3", "BIO19", "BIO18", "PHH2O", "BIO5", "BIO12", "BIO9", "Elevation",
               "BIO1", "BIO17", "BIO13", "BIO4", "BIO2", "SOC", "BIO8", "SRAD", 
               "BIO10", "BIO14", "BIO7", "BIO16", "BIO15", "BIO6", "BIO11"),
  Value = c(9271, 9142, 7029, 4475, 2767, 2688, 2312, 2186,
            1357, 1331, 1234, 1162, 1075, 932, 896, 864, 
            565, 525, 462, 298, 287, 93, 66)
)

# 添加分组信息
data$Group <- "Group 1"
data$Group[data$Variable %in% c("BIO19", "BIO18", "BIO12", "BIO17", "BIO13", "BIO14", "BIO16", "BIO15", "BIO11")] <- "Group 2"
data$Group[data$Variable %in% c("PHH2O", "Elevation", "SOC", "SRAD")] <- "Group 3"

# 绘制条形图
ggplot(data, aes(x = reorder(Variable, Value), y = Value, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Value), hjust = -0.3, size = 3) +  # 在条形图上添加数值标签
  coord_flip() +  # 水平条形图
  theme_bw() +
  scale_fill_manual(values = c("Group 1" = "darkorange", "Group 2" = "steelblue", "Group 3" = "forestgreen")) +
  labs(title = "Variable Values",
       x = "Variable",
       y = "Value",
       fill = "Group") +
  theme(plot.title = element_text(hjust = 0.5))  # 标题居中
