# 加载ggplot2包
library(ggplot2)
library(tidyr)

df = read.csv("191sample_reseq.csv", header = T, row.names = 1)

# 定义列名列表，排除非数值型列'Sample'
# 确认列名列表
columns <- c("Raw_reads", "Raw_bases", "Clean_reads", "Clean_bases", "Effective_rate", "Read_length", "Q20", "Q30", "GC", "Mapped", "Average_sequencing_depth", "Coverage", "Coverage_at_least_4X", "Coverage_at_least_10X", "Coverage_at_least_20X", "SNP_num", "Transition", "Transversion", "Ts_Tv", "Heter_num", "Heter_ratio", "Hom_num", "Hom_ratio")

library(ggplot2)
library(tidyr)

# 将数据转换为长格式
df_long <- pivot_longer(df, cols = columns, names_to = "Variable", values_to = "Value")

# 使用长格式数据绘制图形，每个变量一个图
plots_list <- lapply(columns, function(col_name) {
  ggplot(df_long[df_long$Variable == col_name, ], aes(x = Variable, y = Value)) +
    geom_violin(trim = FALSE) +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), size = 1, alpha = 0.5) +
    labs(x = col_name, y = "") +
    theme_bw() +
    theme(axis.ticks.x = element_blank(), 
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(col_name)
})

# 保存图形或者打印出来
# 这里是打印的例子，你也可以将它们保存到文件中
pdf("output_violin_plots.pdf", width = 4, height = 6)
for (plot in plots_list) {
  print(plot)
}
dev.off()
