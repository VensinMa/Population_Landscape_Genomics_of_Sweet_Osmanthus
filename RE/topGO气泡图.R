
library(dplyr)
library(readxl)

file1 <- read.csv("1861_geneid.csv")
file2 <- read_excel("Osmanthus.genomic.pep.Nr.annotations.xlsx")

# 基于GeneID合并文件，并保留所有文件一中的列和文件二中匹配的行
matched_rows <- left_join(file1, file2, by = "GeneID")

# 查看合并后的前几行数据确保正确合并
head(matched_rows)

# 将matched_rows保存为制表符分隔的文件
write.table(matched_rows, "matched_rows.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)

###############################################################################
library(ggplot2)
library(dplyr)
library(readxl)

# 读取数据
go_data <- read_excel("go富集.xlsx", sheet = "140")

# 过滤出 Pvalue 小于 0.05 的行，然后选择每个 GO 分类的前 10 项
top_terms_per_category <- go_data %>%
  filter(Pvalue < 0.05) %>%
  group_by(GO) %>%
  slice_max(order_by = Significant, n = 10) %>%
  ungroup()
top_terms_per_category <- go_data %>%
  filter(Pvalue < 0.05) %>%
  group_by(GO) %>%
  slice_max(order_by = Significant, n = 10, with_ties = FALSE) %>%
  ungroup()



ggplot(top_terms_per_category, aes(x = reorder(Term, Significant), y = Significant, size = Significant, color = Pvalue)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(3, 6)) +
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~GO, scales = "free_y", ncol = 1) +  # ncol = 1 使得分面垂直排列
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(hjust = 1),
    strip.text.x = element_text(size = 8),
    strip.background = element_blank(),  # 移除分面标签的背景
    panel.grid.major = element_blank(),  # 移除主要网格线
    panel.grid.minor = element_blank(),  # 移除次要网格线
    panel.border = element_blank()       # 移除面板边框
  ) +
  theme_bw() +
  labs(
    x = "",
    y = "GO Term",
    size = "Count",
    color = "P-value"
  ) +
  coord_flip()  # 交换 x 和 y 轴以更好地显示长标签



ggplot(top_terms_per_category, aes(x = Significant, y = reorder(Term, Significant), size = Significant, color = Pvalue)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(3, 6)) +
  scale_color_gradient(low = "blue", high = "red") +
  facet_wrap(~GO, scales = "free_y", ncol = 1) +  # ncol = 1 使得分面垂直排列
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.x = element_text(size = 8),
    #strip.background = element_blank(),  # 移除分面标签的背景
    #panel.grid.major = element_blank(),  # 移除主要网格线
   # panel.grid.minor = element_blank(),  # 移除次要网格线
    #panel.border = element_blank()       # 移除面板边框
  ) +
  labs(
    x = "Significant",
    y = "GO Term",
    size = "Count",
    color = "P-value"
  )


ggplot(top_terms_per_category, aes(x = Significant, y = reorder(Term, Significant), size = Significant, color = Pvalue)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(3, 6)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") + # 使用facet_grid将标签放在侧边
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(hjust = 1),
    strip.text.y = element_text(angle = 0),  # 将分面标签设置为竖直显示
    #strip.background = element_blank(),  # 移除分面标签的背景
    #panel.grid.major = element_blank(),  # 移除主要网格线
    #panel.grid.minor = element_blank(),  # 移除次要网格线
    #panel.border = element_blank(),       # 移除面板边框
    strip.placement = "outside"           # 将分面标签放在外侧
  ) +
  labs(
    x = "Significant",
    y = "GO Term",
    size = "Count",
    color = "P-value"
  )

ggplot(top_terms_per_category, aes(x = Significant, y = reorder(Term, Significant), color = Pvalue, size = Significant)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(3, 6)) +
  scale_color_gradient(low = "red", high = "blue") +
  facet_grid(GO ~ ., scales = "free_y", space = "free_y") + # 使用facet_grid将标签放在侧边
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.x = element_text(angle = 90, hjust = 1, color = "black"),  # X轴文字颜色
    axis.text.y = element_text(hjust = 1, color = "black"),              # Y轴文字颜色
    strip.text.y = element_text(angle = 0, color = "black"),             # 分面标签文字颜色
    strip.background = element_rect(fill = "#D9D9D9", color = "black"),    # 分面标签背景和边框颜色
    panel.border = element_rect(color = "black", fill = NA),             # 面板边框颜色
    legend.text = element_text(color = "black"),                         # 图例文字颜色
    legend.title = element_text(color = "black"),                        # 图例标题颜色
    strip.placement = "outside"                                          # 将分面标签放在外侧
  ) +
  labs(
    x = "Significant",
    y = "GO Term",
    color = "P-value",
    size = "Count"
  )


