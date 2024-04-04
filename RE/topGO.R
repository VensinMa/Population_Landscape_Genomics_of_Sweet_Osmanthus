
# 加载包
rm(list=ls())
library(topGO)
library(Rgraphviz)


getwd()
dir_Manhattan <- "topGO"
# 检查目录是否存在
if (!file.exists(dir_Manhattan)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_Manhattan)
}
setwd(dir_Manhattan)
getwd()

# 设置输入文件
input="140_wxz_uniq.txt"  #待分析基因名称的列表
mapfile="all_gene_GOs.txt"    #所有基因GO注释结果 两列 空格" "作为分隔符

# 开始分析
gene_id = readMappings(file = mapfile,sep = " ", IDsep = ",") #如果是读取其他文件格式，后面参数还需要修改
gene_names = names(gene_id)
my_genes = read.table(input)[,1]

gene_list = rep(1,length(gene_id))
names(gene_list) = names(gene_id)

gene_list[match(my_genes,names(gene_list))] = 0
top_diff_genes = function(allScore){return(allScore<0.01)}


# 1. BP 富集分析
## new() 创建一个 topGO 的对象，然后对这个对象做检验
bp_go = new("topGOdata",
            nodeSize = 6,
            ontology="BP",
            allGenes = gene_list,
            annot = annFUN.gene2GO,
            gene2GO = gene_id,
            geneSel=top_diff_genes)

## 做显著性检验
result_BP.fisher = runTest(bp_go,
                         algorithm = "weight01",
                         statistic = "fisher")

## 提取基因 table
bp_allres = GenTable(bp_go,
                  p.value = result_BP.fisher,
                  orderBy = "p.value",
                  ranksOf = "classic",
                  topNodes = attributes(result_BP.fisher)$geneData[4])

## 生成文件，后面画图都可以用这个表
write.table(bp_allres,
            file = paste(input,".GO_BP_140.xls",sep=""),
            sep="\t", col.names=TRUE, row.names=FALSE)

## 输出矢量图
pdf(paste(input,".BP.pdf",sep=""))
showSigOfNodes(bp_go,
               score(result_BP.fisher),
               firstSigNodes = 10,
               useInfo = "all") #设置节点数量，10个或者20个更多都可以
dev.off()

## 输出像素图
png(paste(input,".BP.png",sep=""))
showSigOfNodes(bp_go, score(result_BP.fisher), firstSigNodes = 10, useInfo = "all")
dev.off()

# 2. MF 富集分析（同理）
## 创建一个 topGO 的对象
mf_go = new("topGOdata",
            nodeSize = 6,
            ontology="MF",
            allGenes = gene_list,
            annot = annFUN.gene2GO,
            gene2GO = gene_id,
            geneSel=top_diff_genes)

## 做检验
result_MF.fisher = runTest(mf_go,
                         algorithm = "weight01",
                         statistic = "fisher")

## 提取基因 table
mf_allres = GenTable(mf_go ,
                  p.value = result_MF.fisher,
                  orderBy = "p.value",
                  ranksOf = "classic",
                  topNodes = attributes(result_MF.fisher)$geneData[4])

## 生成文件，后面画图都可以用这个表
write.table(mf_allres,
            file = paste(input,".GO_MF_140.xls",sep=""),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

## 输出矢量图
pdf(paste(input,".MF.pdf",sep=""))
showSigOfNodes(mf_go,
               score(result_MF.fisher),
               firstSigNodes = 10,
               useInfo = "all") #设置节点数量，10个或者20个更多都可以
dev.off()

## 输出像素图
png(paste(input,".MF.png",sep=""))
showSigOfNodes(mf_go, score(result_KS.elim), firstSigNodes = 10, useInfo = "all")
dev.off()


# 3. CC节点的富集分析(同理)
## 创建一个 topGO 的对象
cc_go = new("topGOdata",
            nodeSize = 6,
            ontology="CC",
            allGenes = gene_list,
            annot = annFUN.gene2GO,
            gene2GO = gene_id,
            geneSel=top_diff_genes)

## 做检验
result_CC.fisher = runTest(cc_go,
                         algorithm = "weight01",
                         statistic = "fisher")

## 提取基因 table
cc_allres = GenTable(cc_go,
                  p.value = result_CC.fisher,
                  orderBy = "p.value",
                  ranksOf = "classic",
                  topNodes = attributes(result_CC.fisher)$geneData[4])

## 生成文件，后面画图都可以用这个表
write.table(cc_allres,
            file = paste(input,".GO_CC_140.xls",sep=""),
            sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

## 输出矢量图
pdf(paste(input,".CC.pdf",sep=""))
showSigOfNodes(cc_go,
               score(result_CC.fisher),
               firstSigNodes = 10,
               useInfo = "all") #设置节点数量，10个或者20个更多都可以
dev.off()

## 输出像素图
png(paste(input,".CC.png",sep=""))
showSigOfNodes(cc_go, score(result_CC.fisher), firstSigNodes = 10, useInfo = "all")
dev.off()


df <- read.csv("140_go.csv", stringsAsFactors = FALSE)
# 使用Benjamini-Hochberg方法校正p值
df$p.adjusted <- p.adjust(df$p.value, method = "BH")
head(df)
# 将校正后的数据写入新的CSV文件
write.csv(df, "140_go_adjusted_p_values.csv", row.names = FALSE)



