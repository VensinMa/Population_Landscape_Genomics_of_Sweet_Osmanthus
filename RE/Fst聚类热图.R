
library(ComplexHeatmap)
library(circlize) 
library(pheatmap)

# 读取数据
fst_matrix <- read.csv("186samples_32pops_fst_matrix.csv", 
                       header = TRUE, row.names = 1)



mycol1<-colorRampPalette(c("navy", "white", "orange"))(100)

# 绘制聚类热图
pheatmap(fst_matrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         scale = "none",
         color = mycol1,
         display_numbers = T,
         fontsize = 8,
         fontsize_row = 6,
         fontsize_col = 7,
         number_color = "black",
         cellwidth = 15,
         cellheight = 8,
         border_color = "white")

