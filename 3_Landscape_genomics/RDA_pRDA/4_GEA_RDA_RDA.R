library(psych)    
library(vegan) 
library(adegenet)
library(LEA)
library(data.table)
library(parallel)
library(cols4all)

###########################  设置工作目录  ####################################
setwd("C:/Rstudio/RStudio/Workspace/GEA_2024")
getwd()
# 目标目录
dir_name <- "RDA"
# 检查目录是否存在
if (!file.exists(dir_name)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_name)
}
# 改变工作目录到目标目录
setwd(dir_name)
getwd()
## [1] "C:/Rstudio/RStudio/Workspace/GEA_2024/RDA"

#可直接加载已插补好的基因型数据  
# load("gen.imp.RData")

# 读取文件
rawdata <- fread("C:/Users/mawenxin/Downloads/Linux_down/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink.raw")
dim(rawdata)   # [1]     194 1497185

# 提取第7列到最后一列  即基因型
rawgeno <- rawdata[, 7:ncol(rawdata)]
dim(rawgeno)   # [1]     194 1497179

# 读取 SNP 标识符文件
snp_ids <- fread("C:/Users/mawenxin/Downloads/Linux_down/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink.snp.id", 
                 header = FALSE)

# 检查 SNP 标识符的数量是否与基因型数据的列数一致
if (ncol(rawgeno) != nrow(snp_ids)) {
  stop("SNP 标识符的数量与基因型数据的列数不一致。")
}

# 将 SNP 标识符文件中的值作为基因型数据的列名
colnames(rawgeno) <- snp_ids$V1

# 检查是否成功替换列名
head(colnames(rawgeno))


# 将提取后的数据保存到新文件
fwrite(rawgeno, "194samples_filtered.LD.pruned.geno.raw", 
       sep = "\t", quote = FALSE)

## 基因型插补
gen.imp <- apply(rawgeno, 2, function(x) {
  if(any(is.na(x))) {
    # 计算最常见的基因型
    modes <- as.numeric(names(which.max(table(x))))
    # 插补缺失值
    x[is.na(x)] <- modes
  }
  return(x)
})

save(gen.imp, file = "gen.imp.RData")
# load("gen.imp.RData")


## 加载群体经纬度环境数据
env_data = read.csv("C:/Rstudio/RStudio/Workspace/GEA_2024/extracted_data/194sample_id_group_lat_lon_ele_with_env_vars.csv",
                    header = T, row.names = 1)

pairs.panels(env_data[, 4:26], smooth = TRUE, scale = TRUE, density = TRUE,
             method = "pearson", pch = 20, lm = FALSE, cor = TRUE,
             jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE,
             rug = TRUE, breaks = "Sturges",  wt = NULL, ellipses = TRUE, 
             cex.cor = 5, cex = 1, smoother = FALSE, stars = TRUE, alpha = .05, 
             hist.border = "black")

pairs.panels(env_data[, c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                          "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")], 
             smooth = TRUE, scale = TRUE, density = TRUE,
             method = "pearson", pch = 20, lm = FALSE, cor = TRUE,
             jiggle = FALSE, factor = 2, hist.col = "cyan", show.points = TRUE,
             rug = TRUE, breaks = "Sturges",  wt = NULL, ellipses = TRUE, 
             cex.cor = 5, cex = 1, smoother = FALSE, stars = TRUE, alpha = .05, 
             hist.border = "black")

## pdf导出尺寸为16 * 16

# load("of.rda.RData")
# 假设pred是你的环境变量数据框，gen.imp是响应变量矩阵
of.rda <- rda(gen.imp ~ BIO2 + BIO8 + BIO9 + BIO10 + BIO12 + BIO15 + 
                BIO17 + BIO18 + SRAD + SOC + PHH2O, 
              env_data[, c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                           "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")])
# 保存多个对象到一个文件
# save(of.rda, file = "of.rda.RData")

# 查看RDA结果
summary(of.rda)
screeplot(of.rda)
vif.cca(of.rda)
RsquareAdj(of.rda)
A = summary(eigenvals(of.rda, model = "constrained"))
write.csv(A, file = "Importance_of_components_eigenvals_RDA.csv")

png("RDA_AllLoci.png", width=15, height=15, units='in', res=500)
plot(of.rda, scaling = 1,  # 'n'表示不绘制任何东西，仅设置绘图区域
     cex.axis = 1.4, 
     cex.lab = 1.4, 
     xlab = "RDA1 site score", 
     ylab = "RDA2 site score", 
     las = 1)  # 设置轴标签方向
dev.off()
??plot

c4a_gui()
library(cols4all)
set.seed(11)
colors = c4a("dynamic", length(unique(env_data$group)))

# 假定colors包含了三个颜色代码
colors <- c("#F47F72", "#FFD700", "#7FB2D5")

# 设置env_data$pop的水平
levels(env_data$pop) <- c("Central", "East", "West")

# 创建一个命名向量，将分组名称映射到颜色
color_mapping <- setNames(colors, levels(env_data$group))

# 将env_data$group中的每个分组名称替换为对应的颜色
env_data$group_color <- color_mapping[env_data$group]

# 打印结果，查看分组对应的颜色
bg = env_data$group_color


# bg <- c("#ff7f00", "#1f78b4", "#ffff33", "#a6cee3", "#33a02c")

png("RDA_AllLoci_final.png", width=10, height=10, units='in', res=500)
# axes 1 & 2
plot(of.rda, type="n", scaling=1)
#points(of.rda, display="species", pch=20, cex=0.7, col="gray32", scaling=1)           # the SNPs
points(of.rda, display="sites", pch=21, cex=1.3, col="gray20", scaling=1, bg=bg) # the samples
text(of.rda, scaling=1, display="bp", col="#0868ac", cex=1)                           # the predictors
legend("bottom", legend=levels(eco), bty="n", col="black", pch=21, cex=1, pt.bg=colors)
dev.off()
??legend



load.rda <- scores(of.rda, choices=c(1:3), display="species")  # Species scores for the first three constrained axes

hist(load.rda[,1], main="Loadings on RDA1")
hist(load.rda[,2], main="Loadings on RDA2")
hist(load.rda[,3], main="Loadings on RDA3") 
#hist(load.rda[,4], main="Loadings on RDA4")
#hist(load.rda[,5], main="Loadings on RDA5")
#hist(load.rda[,6], main="Loadings on RDA6") 

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}

cand1 <- outliers(load.rda[,1],3) # 21890   30779
cand2 <- outliers(load.rda[,2],3) # 16000   22958
cand3 <- outliers(load.rda[,3],3) # 15405   19562
length(cand1)
length(cand2)
length(cand3)
names_list <- list(names(cand1), names(cand2), names(cand3))

# 使用Reduce函数配合union一次性合并所有names，自动去除重复项
combined_cand <- Reduce(union, names_list)
ncand = length(combined_cand)
ncand  # [1] 51243

# 写入CSV文件
write.csv(data.frame(names = combined_cand), "RDA_cand.csv", row.names = FALSE)



cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)
length(cand$snp) # 53295
# 去除 cand 中的重复项
unique_cand <- unique(cand$snp)
length(unique_cand)  # 这应该与 ncand 相同


 foo <- matrix(nrow=ncand, ncol=23)  # 23 columns for 23 environmental variables
colnames(foo) <- c("BIO1", "BIO2", "BIO3", "BIO4", "BIO5", "BIO6", "BIO7", "BIO8", "BIO9", "BIO10", 
                   "BIO11", "BIO12", "BIO13", "BIO14", "BIO15", "BIO16", "BIO17", "BIO18", "BIO19", 
                   "Elevation", "SRAD", "SOC", "PHH2O")


for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- gen.imp[,nam]
  foo[i,] <- apply(env_data[, 4:26],2,function(x) cor(x,snp.gen))
}

cand <- cbind.data.frame(cand,foo)  
head(cand)

length(cand$snp[duplicated(cand$snp)])  # 7 duplicate detections
foo <- cbind(cand$axis, duplicated(cand$snp)) 
table(foo[foo[,1]==1,2]) # no duplicates on axis 1
table(foo[foo[,1]==2,2]) #  7 duplicates on axis 2
table(foo[foo[,1]==3,2]) # no duplicates on axis 3

cand <- cand[!duplicated(cand$snp),] # remove duplicate detections


for (i in 1:length(cand$snp)) {
  bar <- cand[i,]
  cand[i,12] <- names(which.max(abs(bar[4:11]))) # gives the variable
  cand[i,13] <- max(abs(bar[4:11]))              # gives the correlation
}

colnames(cand)[12] <- "predictor"
colnames(cand)[13] <- "correlation"

table(cand$predictor) 

