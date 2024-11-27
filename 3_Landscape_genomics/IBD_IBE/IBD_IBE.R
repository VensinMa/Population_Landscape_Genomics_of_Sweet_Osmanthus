# 加载所需的R包
library(vegan)
library(geosphere)
library(ggplot2)
library(cowplot)

getwd()
# 设置工作目录
setwd("C:/RStudio/RStudio/Workspace/IBD_IBE")

# 读取环境数据和Fst数据
env <- read.csv("input/32pop_means_env_vars.csv", head = TRUE, row.names = 1)  # 读取环境数据
fst <- read.csv("input/179w_fst_matrix.csv", head = TRUE, row.names = 1)  # 读取Fst数据
# 确保env中包含fst所有的行
if(all(row.names(fst) %in% row.names(env))) {
  # 使用fst的行名来重新排序env
  env_reordered <- env[row.names(fst), ]
} else {
  # 如果env中不包含fst所有的行，打印警告信息
  print("Warning: Not all rows in fst are present in env.")
}

# 检查重新排序后的env
head(env_reordered)
env = env_reordered

# 提取环境数据的子集并进行标准化
ENV <- as.data.frame(env[, c(3:ncol(env))])
ENV <- scale(ENV, center = TRUE, scale = TRUE)
rownames(ENV) <- row.names(env)
rownames(ENV)

# 计算环境距离
env_dist <- vegdist(ENV, method = "euclidean", binary = FALSE, diag = FALSE, 
                    upper = FALSE, na.rm = FALSE)
write.csv(as.matrix(env_dist), file = "env_dist.csv", quote = FALSE)  # 将环境距离保存为CSV文件

# 提取地理距离数据
dist <- as.data.frame(env[, c("lon","lat")])
rownames(dist) <- row.names(env)
rownames(dist)

# 计算地理距离
muer.dists <- distm(dist, fun = distVincentyEllipsoid)
rownames(muer.dists) <- row.names(env)
colnames(muer.dists) <- row.names(env)
write.csv(as.matrix(muer.dists), file = "geo_dist.csv", quote = FALSE)  # 将地理距离保存为CSV文件

# 计算环境距离与地理距离之间的相关性
env_dist <- read.csv("env_dist.csv", head = TRUE, row.names = 1)
mantel(muer.dists, env_dist, method = "pearson", permutations = 999)  # 全部数据

# 计算Fst和地理距离之间的相关性  IBD
mantel(fst, muer.dists, method = "pearson", permutations = 999)  # 全部数据

# 计算Fst和环境距离之间的相关性  IBE
#fst <- read.csv("input/1861_fst_matrix.csv", head = TRUE, row.names = 1)
mantel(fst, env_dist, method = "pearson", permutations = 999)  # 全部数据

# 计算Fst、环境距离和地理距离之间的相关性（部分相关性） IBE partial
mantel.partial(fst, env_dist, muer.dists, method = "pearson", permutations = 999)  # 全部数据

# 计算Fst、地理距离和环境距离之间的相关性（部分相关性） IBD partial
mantel.partial(fst,muer.dists,env_dist,method="pearson",permutations=999)

##################################  PLOT  ######################################
library(ggplot2)
library(cowplot)

geo_dit=read.csv("geo_dist.csv",head=F)
env_dist=read.csv("env_dist.csv",head=F)
FST_adaption=read.csv("input/1861_fst_matrix.csv",head=F)
FST_Neutral=read.csv("input/179w_fst_matrix.csv",head=F)
trans <- function(raw_data){
  raw_data=raw_data[-1,-1]
  out_data=data.frame(raw_data[,1])
  colnames(out_data)="value"
  for (i in 2:32){
    temp=data.frame(raw_data[i:32,i])
    colnames(temp)="value"
    out_data=rbind(out_data,temp)
  }
  out_data=na.omit(out_data)
  return(out_data)
}

plot_data <- data.frame(
  geo_dit = trans(geo_dit)$value,
  env_dist = trans(env_dist)$value,
  fst1779 = trans(FST_adaption)$value,
  fst53w = trans(FST_Neutral)$value
)

colnames(plot_data)=c("geo_dit","env_dist","FST_adaption","FST_Neutral")
write.csv(plot_data,file="plot_data.csv",quote=F,row.names = F)

plot_data=read.csv("plot_data.csv",head=T)
p1=ggplot(plot_data)+
  geom_point(aes(x=(geo_dit/1000),y=FST_adaption),size = 3,alpha=0.9,color="gray20",shape=21,fill="#F47F72")+
  geom_point(aes(x=(geo_dit/1000),y=FST_Neutral),size = 3,alpha=0.9,color="gray20",shape=21,fill="#7FB2D5")+
  geom_smooth(aes(x=(geo_dit/1000), y=FST_Neutral),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="gray70", fill="#b8c0c8",linewidth = 1.5,fullrange = F) +
  geom_smooth(aes(x=(geo_dit/1000), y=FST_adaption),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#E8ACB7",fill="#b8c0c8",fullrange = F, linewidth = 1.5) +
  labs(x = "Geographical Distance (100km)",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.6, linetype = 1), panel.background = element_blank())+
  theme_bw()+
  theme(text=element_text(family="sans"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))

p1
p2=ggplot(plot_data)+
  geom_point(aes(x=env_dist,y=FST_adaption),size = 3,alpha=0.9,color="gray20",shape=21,fill="#F47F72")+
  geom_point(aes(x=env_dist,y=FST_Neutral),size = 3,alpha=0.9,color="gray20",shape=21,fill="#7FB2D5")+
  geom_smooth(aes(x=env_dist, y=FST_Neutral),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="gray70", fill="#b8c0c8",linewidth = 1.5,fullrange = F) +
  geom_smooth(aes(x=env_dist, y=FST_adaption),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#E8ACB7",fill="#b8c0c8",fullrange = F, linewidth = 1.5) +
  labs(x = "Environment Distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5)+
  theme(panel.border = element_rect(colour = "black", linewidth = 0.6, linetype = 1), panel.background = element_blank())+
  theme_bw()+
  theme(text=element_text(family="sans"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        axis.line.x=element_line(colour="black"),
        axis.line.y=element_line(colour="black"),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))
p2
all=plot_grid(p1,p2,align ="v",labels=c("c","d"),label_size = 20,label_fontfamily = "sans",label_fontface = 1,ncol=1)
ggsave(all,file="P32.pdf",width=5.5,height=8)

all
