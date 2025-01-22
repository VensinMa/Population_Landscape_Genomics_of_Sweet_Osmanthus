library(psych)    
library(vegan) 
library(adegenet)
library(LEA)
library(data.table)
library(parallel)
library(sp) # 用于空间点数据处理
library(caret) # install.packages("caret", dependencies = c("Depends", "Suggests"))
library(lavaan)
library(ggpubr)
library(adespatial)
#source ('https://raw.githubusercontent.com/zdealveindy/anadat-r/master/scripts/NumEcolR2/quickMEM.R')

###########################  设置工作目录  ####################################
getwd()
setwd("/public1/guop/mawx/gc/R/GEA")
getwd()
#可直接加载已插补好的基因型数据  
# load("gen.imp.RData")

###  直接读取 lfmm格式文件
rawgeno <- fread("input/311_filtered.recode_noscaffold.plink.recodeA.lfmm")

# 读取 SNP 标识符文件
snp_ids <- fread("input/311_filtered.recode_noscaffold.plink.ID", 
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
#fwrite(rawgeno, "194samples.filtered.missing0.05.raw", sep = "\t", quote = FALSE)

## 基因型插补  lfmm文件缺失基因型为9
gen.imp <- apply(rawgeno, 2, function(x) {
  # 检查是否有 9（表示缺失值）
  if (any(x == 9)) {
    # 计算列中不为9的最常见基因型
    modes <- as.numeric(names(which.max(table(x[x != 9]))))
    # 替换为常见基因型
    x[x == 9] <- modes
  }
  return(x)
})

save(gen.imp, file = "gen.imp.3284w.RData")
#load("gen.imp.RData")

## 加载群体经纬度环境数据
env_data = read.csv("input/Climate_current_311samples.csv",
                    header = T, row.names = 1)
new_order = read.table("input/311_filtered.recode_noscaffold.plink.recodeA.nosex", header = F)
new_order = new_order$V1
env_data = env_data[new_order, ]
str(env_data)

# load("of.rda.RData")
# 假设pred是你的环境变量数据框，gen.imp是响应变量矩阵
of.rda <- rda(gen.imp ~ BIO2 + BIO3 + BIO4 + BIO8 + BIO10 + BIO12 + 
                BIO15 + BIO18 + SRAD + SOC + Elevation, 
              env_data[, c("BIO2", "BIO3", "BIO4", "BIO8", "BIO10", "BIO12", 
                           "BIO15", "BIO18", "SRAD", "SOC", "Elevation")])
# 保存of.rda
save(of.rda, file = "of.rda.3284w.RData")

# 查看RDA结果
screeplot(of.rda)
summary(of.rda)
vif.cca(of.rda)
RsquareAdj(of.rda)
summary(eigenvals(of.rda, model = "constrained"))

sink("Result_RDA.txt")
summary(of.rda)
vif.cca(of.rda)
RsquareAdj(of.rda)
summary(eigenvals(of.rda, model = "constrained"))
sink()

write.csv(summary(eigenvals(of.rda, model = "constrained")), file = "Importance_of_components_eigenvals_RDA.csv")

## png("RDA_AllLoci.png", width=15, height=15, units='in', res=500)
## plot(of.rda, scaling = 1,  # 'n'表示不绘制任何东西，仅设置绘图区域
##     cex.axis = 1.4, 
##     cex.lab = 1.4, 
##     xlab = "RDA1 site score", 
##    ylab = "RDA2 site score", 
##     las = 1)  # 设置轴标签方向
##dev.off()

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

cand1 <- outliers(load.rda[,1],3) # 26468
cand2 <- outliers(load.rda[,2],3) # 31786
cand3 <- outliers(load.rda[,3],3) # 31053
length(cand1)
length(cand2)
length(cand3)
names_list <- list(names(cand1), names(cand2), names(cand3))

# 使用Reduce函数配合union一次性合并所有names，自动去除重复项
combined_cand <- Reduce(union, names_list)
ncand = length(combined_cand)
ncand  # [1] 89262

# 写入CSV文件
write.csv(data.frame(names = combined_cand), "RDA_cand.csv", row.names = FALSE)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))

colnames(cand1) <- colnames(cand2) <- colnames(cand3) <- c("axis","snp","loading")

cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)
length(cand$snp) # 89307
unique_cand <- unique(cand$snp) # 去除 cand 中的重复项
length(unique_cand)  # 这应该与 ncand 相同  # 89262


###################################           #################################

env = read.csv("input/Climate_current_311samples.csv",
               header = T, row.names = 1)

env_data = env[, c("BIO2", "BIO3", "BIO4", "BIO8", "BIO10", "BIO12", 
                   "BIO15", "BIO18", "SRAD", "SOC", "Elevation")]
env_select = as.matrix(env_data)
head(env_select)

# 生成dbMEM
dbmem_data <- dbmem(env[, c("lon", "lat")])
dbmem = as.matrix(dbmem_data)
dbmen.rda <- rda(env_select ~ ., dbmem_data)
sink("RDA.txt")
anova.cca(dbmen.rda, permutations = 999)
## 模型显著 可进行正向选择
dbmen.rda.r2a <- RsquareAdj(dbmen.rda)$adj.r.squared
dbmen.rda.r2a # 0.898679

## 正向选择
dbmen.fwd <- forward.sel(env_select, dbmem, adjR2thresh = dbmen.rda.r2a)
nrow(dbmen.fwd) # 22

## 排序
dbmen.fwd.sign.sort <- sort(dbmen.fwd$order)
dbmen.fwd.sign.sort

## 提取显著dbmen
dbmen.fwd.sign <- dbmem[, dbmen.fwd.sign.sort]

################################      RDA    ###################################
of_rda <- rda(gen.imp ~ env_select)

of_partial_rda <- rda(gen.imp ~ env_select + Condition(dbmen.fwd.sign))
summary(of_partial_rda)

of_partial_rda_r2 <- RsquareAdj(of_partial_rda)$adj.r.squared
of_partial_rda_r2
#save(of_partial_rda, file = "1861_of_partial_rda.RData")
# load("1861_of_partial_rda.RData")

vp <- varpart(gen.imp, env_select, dbmen.fwd.sign)
#save(vp, file = "1861_of_partial_rda_vp.RData")
# load("1861_of_partial_rda_vp.RData")
print(vp)

sink("pvalue_RDA_pRDA.607w.txt")
# 测试 F ~ env 的显著性
rda_env <- rda(gen.imp ~ env_select)
anova(rda_env, permutations = 99)
# 测试 F ~ geo 的显著性
rda_geo <- rda(gen.imp ~ dbmen.fwd.sign)
anova(rda_geo, permutations = 99)
# 测试 F ~ env + geo 的显著性    env + geo的组合贡献
rda_env_geo <- rda(gen.imp ~ env_select + dbmen.fwd.sign)
anova(rda_env_geo, permutations = 99)

# 测试 F ~ env | geo 的显著性
prda_env <- rda(gen.imp ~ env_select + Condition(dbmen.fwd.sign))
anova(prda_env, permutations = 99)
# 测试 F ~ geo | env 的显著性
prda_geo <- rda(gen.imp ~ dbmen.fwd.sign + Condition(env_select))
anova(prda_geo, permutations = 99)
print(vp)

sink()

rda_model_env <- rda(gen.imp ~ BIO2 + BIO8 + BIO9 + BIO10 + BIO12 + BIO15 + 
                       BIO17 + BIO18 + SRAD + SOC + PHH2O + Condition(dbmen.fwd.sign),
                     data = env_data)
anova_res_env <- anova.cca(rda_model_env, permutations = 999)
print(anova_res_env)
anova_res_env <- anova.cca(rda_model_env, permutations = 999,by = "margin")
print(anova_res_env)

rda_model_spatial <- rda(gen.imp ~ dbmen.fwd.sign)
anova_res_spatial <- anova.cca(rda_model_spatial, permutations = 999)
print(anova_res_spatial)

rda_model_combined <- rda(gen.imp ~ env_select + dbmen.fwd.sign)
anova_res_combined <- anova.cca(rda_model_combined, permutations = 999)
print(anova_res_combined)

rda_model_env_partial <- rda(gen.imp ~ BIO2 + BIO8 + BIO9 + BIO10 + BIO12 + BIO15 + 
                               BIO17 + BIO18 + SRAD + SOC + PHH2O + Condition(dbmen.fwd.sign),
                             data = env_data)
anova_res_env_partial <- anova.cca(rda_model_env_partial, permutations = 999)
print(anova_res_env_partial)

rda_model_spatial_partial <- rda(gen.imp ~ dbmem + Condition(env_select))
anova_res_spatial_partial <- anova.cca(rda_model_spatial_partial, permutations = 999)
print(anova_res_spatial_partial)
# 使用999次排列来评估模型中每个变量的显著性
anova_res_terms <- anova.cca(of_partial_rda, by = "axis",permutations = 999,parallel = 16)

# 查看每个变量的显著性结果
print(anova_res_terms)

# 打印每个分量的置换检验结果
print(anova_res_env)
print(anova_res_spatial)



png("venn_partitions.png", width = 6, height = 6, units="in", res=600)
plot(vp)
dev.off()

# 创建一个数据框以保存结果
varpart_results <- data.frame(
  Fraction = c("Environmental", "Spatial", "Environmental|Spatial", "Spatial|Environmental", "Unexplained"),
  Adj.R.squared = c(0.11941, 0.11569, 0.06020, 0.05649, 0.82411),
  Testable = c(TRUE, TRUE, TRUE, TRUE, FALSE)
)

# 打印结果
print(varpart_results)

