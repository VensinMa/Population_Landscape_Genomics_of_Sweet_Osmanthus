getwd()
dir.create("C:/Rstudio/RStudio/Workspace/GEA_2024/gggenes")
setwd("C:/Rstudio/RStudio/Workspace/GEA_2024/gggenes")
library(gggenes)
library(ggplot2)
gggenes_data = read.csv("gggenes_data.csv",header = T)
gggenes_data2 = read.csv("gggenes_data_Superscaffold11.csv",header = T)

ggplot(
  gggenes_data,
  aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  geom_gene_label(align = "left") +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = rep("#FB8072", nrow(gggenes_data)))  +
  theme_genes()

ggplot(
  gggenes_data2,
  aes(xmin = start, xmax = end, y = molecule, fill = gene, label = gene)) +
  geom_gene_arrow(arrowhead_height = unit(6, "mm"), arrowhead_width = unit(4, "mm")) +
  geom_gene_label(align = "left") +
  facet_wrap(~ molecule, scales = "free", ncol = 1) +
  scale_fill_manual(values = rep("#FB8072", nrow(gggenes_data)))  +
  theme_genes()
