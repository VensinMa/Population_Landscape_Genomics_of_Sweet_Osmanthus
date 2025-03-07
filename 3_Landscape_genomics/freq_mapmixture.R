library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(mapmixture)
library(readr)

# 读取数据
freqdata <- read.csv("input/combined_population_freq.csv") %>%
   select(-FREQ)
freqdata_t <- as.data.frame(t(freqdata))  # 转置并转换为数据框
colnames(freqdata_t) <- freqdata_t[1, ]  # 将第一行作为列名
freqdata_t <- freqdata_t[-1, ]  # 删除原来的列名行
freqdata_t <- rownames_to_column(freqdata_t, var = "pop")
coordinates <- read.csv("input/32pop_means_env_vars.csv")[, c("pop", "lat","lon")]
colnames(coordinate) <- c("Site", "Lat", "Lon")
merged_data <- left_join(freqdata_t, coordinate, by = c("pop" = "Site"))
#write.csv(merged_data,file = "output/merged_data.csv")
# 指定要提取的SNP位点列名
snp_columns <- c("Superscaffold3:3313254","Superscaffold5:1278318")
# 提取每个SNP位点的数据并保存到单独的文件中
for (snp_column in snp_columns) {
  site_data <- freqdata_t %>%
    select(pop, all_of(snp_column)) %>%
    mutate(Allele1 = as.numeric(!!sym(snp_column)),  # 确保是数值型
           Allele2 = 1 - as.numeric(!!sym(snp_column))) %>%
    slice(rep(1:n(), each = 5)) %>%
    group_by(pop) %>%
    mutate(Ind = paste0(pop, rep(1:5, each = n() / 5))) %>%
    ungroup() %>%
    select(Site = pop, Ind, Allele1, Allele2)
  # 将 : 替换为 _ 以便用作文件名
  snp_column_filename <- gsub(":", "_", snp_column)
  # 保存到CSV文件，以SNP名称（替换后的）作为文件名
  write_csv(site_data, paste0("output/",snp_column_filename,"_allele_freq",  ".csv"))
}

# 替换1为0.9999，0为0.0001
replace_values <- function(data) {
  data[data == 1] <- 0.9999
  data[data == 0] <- 0.0001
  return(data)
}

# 定义绘图函数
plot_map <- function(data, cluster_cols, cluster_names, plot_title) {
  mapmixture(data, coordinates, crs = 4326,
             cluster_cols = cluster_cols,
             cluster_names = cluster_names,
             boundary = c(xmin = 97, xmax = 123, ymin = 18, ymax = 33),
             pie_size = 0.5,
             land_colour = "#d9d9d9",
             sea_colour = "#deebf7",
             plot_title = plot_title, # 添加标题
             plot_title_size = 14,    # 标题的字体大小
             axis_title_size = 12,    # 轴标题的字体大小
             axis_text_size = 10)     # 轴文本的字体大小
}

## 
Superscaffold3_3313254 = read.csv("output/Superscaffold3_3313254_allele_freq.csv")
Superscaffold5_1278318 = read.csv("output/Superscaffold5_1278318_allele_freq.csv")
Superscaffold11_5332072 = read.csv("output/Superscaffold11_5332072_allele_freq.csv")
## 
c("#D3FFBE","#B4E17A",,"#A4D874","#FEE0A4","#FECB55","#E69800","#F79C42",
  "#E1E1E1","#FFFFFF","#7A9D3D")
plot_map(Superscaffold3_3313254, c("#5AAFCA", "#F16481"), c("Allele1", "Allele2"), "Superscaffold3_3313254")
plot_map(Superscaffold5_1278318, c("#B4E17A", "#F4A556"), c("Allele1", "Allele2"), "Superscaffold5_1278318")

plot_map(Superscaffold11_5332072, c("#6C7DFE", "#FF930B"), c("Allele1", "Allele2"), "Superscaffold11_5332072")

###########################   批量绘图保存    ##################################
# 定义 SNP 位点列表
snp_columns <- c("Superscaffold3:3313254", "Superscaffold3:3313471","Superscaffold5:1278318",
                 "Superscaffold8:24459601","Superscaffold11:5331878","Superscaffold11:5332072",
                 "Superscaffold20:686440")
for (snp_column in snp_columns) {
  site_data <- freqdata_t %>%
    select(pop, all_of(snp_column)) %>%
    mutate(Allele1 = as.numeric(!!sym(snp_column)),  # 确保是数值型
           Allele2 = 1 - as.numeric(!!sym(snp_column))) %>%
    slice(rep(1:n(), each = 5)) %>%
    group_by(pop) %>%
    mutate(Ind = paste0(pop, rep(1:5, each = n() / 5))) %>%
    ungroup() %>%
    select(Site = pop, Ind, Allele1, Allele2)
  # 将 : 替换为 _ 以便用作文件名
  snp_column_filename <- gsub(":", "_", snp_column)
  # 保存到CSV文件，以SNP名称（替换后的）作为文件名
  write_csv(site_data, paste0("output/",snp_column_filename,"_allele_freq",  ".csv"))
}

# 定义 SNP 位点列表
snp_columns <- c("Superscaffold3:3313254", "Superscaffold3:3313471","Superscaffold5:1278318",
                 "Superscaffold8:24459601","Superscaffold11:5331878","Superscaffold11:5332072",
                 "Superscaffold20:686440")

# 替换 ":" 为 "_"
snp_columns <- gsub(":", "_", snp_columns)

# 循环遍历每个 SNP 位点
for (snp_column in snp_columns) {
  # 读取每个 SNP 位点对应的频率数据
  snp_data <- read.csv(paste0("output/", snp_column, "_allele_freq.csv"))
  plot_obj <- plot_map(snp_data, c("#5AAFCA", "#F16481"), c("Allele1", "Allele2"), snp_column)
  ggsave(paste0("picture/", snp_column, "_freq_map.pdf"), plot = plot_obj, width = 7, height = 5)
}
