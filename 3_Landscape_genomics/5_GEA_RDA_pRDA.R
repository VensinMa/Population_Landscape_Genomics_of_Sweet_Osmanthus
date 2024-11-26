library(psych)    
library(vegan) 
library(adegenet)
library(LEA)
library(data.table)
library(parallel)
library(sp) 
library(rgdal) 
library(caret)
library(lavaan)
library(ggpubr)
library(adespatial)


###########################  设置工作目录  ####################################
getwd()
# 目标目录
dir_name <- "RDA&pRDA"
# 检查目录是否存在
if (!file.exists(dir_name)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_name)
}
# 改变工作目录到目标目录
setwd(dir_name)
getwd()

# 读取文件
rawdata <- fread("C:/Users/mawenxin/Downloads/Linux_down/194samples_snp.nounanchor.renamed.filtered.vcftools.LD.pruned.plink.raw")
dim(rawdata)   # [1]     194 1497185

# 提取第7列到最后一列
# 提取第7列到最后一列  即基因型
rawgeno <- rawdata[, 7:ncol(rawdata)]
dim(rawgeno)   # [1]     194 1497179

# 如果需要，可以将提取后的数据保存到新文件
fwrite(rawgeno, "186_filtered.LD.pruned.noContig.geno.raw", 
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
load("gen.imp.RData")
# load("gen.imp.RData")

env = read.csv("186sample_id_POP_lat_lon_ele_with_env_vars.csv",
               header = T, row.names = 1)

env_data = env[, c("BIO2", "BIO8", "BIO9", "BIO10", "BIO12", 
                   "BIO15", "BIO17", "BIO18", "SRAD", "SOC", "PHH2O")]
env_select = as.matrix(env_data)


# 生成dbMEM
dbmem_data <- dbmem(env[, c("lon", "lat")])
dbmem = as.matrix(dbmem_data)
dbmen.rda <- rda(env_select ~ ., dbmem_data)
anova.cca(dbmen.rda, permutations = 999)
## 模型显著 可进行正向选择
dbmen.rda.r2a <- RsquareAdj(dbmen.rda)$adj.r.squared
dbmen.rda.r2a

## 正向选择
dbmen.fwd <- forward.sel(env_select, dbmem, adjR2thresh = dbmen.rda.r2a)
nrow(dbmen.fwd) # 10

## 排序
dbmen.fwd.sign.sort <- sort(dbmen.fwd$order)
dbmen.fwd.sign.sort

## 提取显著dbmen
dbmen.fwd.sign <- dbmem[, dbmen.fwd.sign.sort]

################################      RDA    ###################################

of_partial_rda <- rda(gen.imp ~ BIO2 + BIO8 + BIO9 + BIO10 + BIO12 + BIO15 + 
                        BIO17 + BIO18 + SRAD + SOC + PHH2O + Condition(dbmen.fwd.sign),
                      data = env_data)
summary(of_partial_rda)

save(of_partial_rda, file = "allsnp_of_partial_rda.RData")
#load("of_partial_rda.RData")

vp <- varpart(gen.imp, env_select, dbmen.fwd.sign)
save(vp, file = "of_partial_rda_vp.RData")
load("of_partial_rda_vp.RData")
print(vp)


# 使用999次排列来评估模型中每个变量的显著性
anova_res_terms <- anova.cca(of_partial_rda, by = "terms", permutations = 99, 
                             parallel = 16)

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






