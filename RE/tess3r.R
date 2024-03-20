devtools::install_github("bcm-uga/TESS3_encho_sen")
install.packages("rworldxtra")
# 加载必要的R包
library(tess3r)
library(maps)
library(raster)
library(ggplot2)
library(rworldmap)
library(LEA)
library(data.table)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(rworldxtra)

getwd()
# 目标目录
dir_name <- "tess3r"
# 检查目录是否存在
if (!file.exists(dir_name)) {
  # 如果目录不存在，则创建目录
  dir.create(dir_name)
}
# 改变工作目录到目标目录
setwd(dir_name)
getwd()

# vcf2lfmm("Z:/root/workspace/186sample/186.snpEffAnno.4dtv.LD.vcf", "186_filtered.LD.4DTV.noContig.lfmm", force = TRUE)

# 加载数据
genotype <- fread("Z:/root/workspace/186sample/186.snpEffAnno.4dtv.LD.lfmm")
# 替换基因型矩阵中的 9 为 NA
genotype[genotype == 9] <- NA
coordinates <- read.csv(file = "186sample_lat_lon.csv", header = F,
                        row.names = 1) # 地理坐标
names(coordinates) <- c('V1', 'V2')
coordinates = as.matrix(coordinates)
# 定义裁剪区域的多边形
polyWKT <- "POLYGON((96 33, 124 33, 124 18, 96 18, 96 33))"
clip_poly <- st_as_sfc(polyWKT, crs = st_crs(4326)) # 创建多边形对象，使用WGS 84坐标系统

# 加载世界地图数据
world <- ne_countries(scale = "medium", returnclass = "sf")

# 裁剪地图数据以仅包含指定多边形区域
world_clipped <- st_intersection(world, clip_poly)

# 创建地图
ggplot() +
  geom_sf(data = world_clipped, fill = "lightgreen", color = "black") + # 绘制裁剪后的国家边界和填充颜色
  geom_point(data = data.frame(coordinates), aes(x = coordinates[,1], y = coordinates[,2]), color = "blue", size = 1) + # 在地图上绘制地理坐标点
  labs(x = "Longitude (°E)", y = "Latitude (°N)") + # 添加坐标轴标签
  theme_bw() + # 使用简洁的主题
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # 移除网格线
  xlim(c(96, 124)) + ylim(c(18, 33)) # 限制地图的显示范围以匹配裁剪区域

plot(coordinates, pch = 19, cex = .5, 
     xlab = "Longitude (°E)", ylab = "Latitude (°N)")
map(add = T, interior = F)

# 使用tess3估算祖先系数
# X是基因型矩阵，coord是地理坐标，K是祖先种群数范围，openMP.core.num是使用的CPU核心数
tess3.obj <- tess3(X = genotype, coord = coordinates, K = 1:8, 
                   method = "projected.ls", ploidy = 2, openMP.core.num = 16) 

# 绘制交叉验证分数图，帮助选择K值
plot(tess3.obj, pch = 19, col = "blue",
     xlab = "Number of ancestral populations",
     ylab = "Cross-validation score")

# 为K=5个种群提取Q矩阵并绘制条形图表示祖先比例
q.matrix <- qmatrix(tess3.obj, K = 5)
q.matrix <- read.table(file = "186_filtered.LD.pruned.4DTV.noContig.5.Q",sep = " ")
q.matrix <- as.matrix(q.matrix)
class(q.matrix) <- 'tess3Q'


# 绘制条形图并捕获返回值
bp <- barplot(q.matrix, border = NA, space = 0, 
              xlab = "Individuals", ylab = "Ancestry proportions", 
              main = "Ancestry matrix")

# 添加自定义的 X 轴标签
axis(1, at = 1:nrow(q.matrix), labels = bp$order, las = 3, cex.axis = .4) 
# 自定义颜色调色板
my.colors <- c("#66000D", "#00441B", "#08306B", "wheat","olivedrab")
my.palette <- CreatePalette(my.colors, 9)
# 使用自定义颜色再次绘制条形图
barplot(q.matrix, border = NA, space = 0, 
        main = "Ancestry matrix", 
        xlab = "Individuals", ylab = "Ancestry proportions", 
        col.palette = my.palette)

# 在地理地图上插值并绘制Q矩阵的值
raster_data <- ("D:/DATA/worldclim/Historical climate data/wc2.1_2.5m_elev/wc2.1_2.5m_elev_cn.tif")  
plot(q.matrix, coordinates, method = "map.max", 
     interpol = FieldsKrigModel(10),window = c(95, 125,20, 35),
     raster.filename = raster_data,
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), col.palette = my.palette)

ggtess3Q(q.matrix, coordinates, resolution = c(300, 300), window = c(95, 125,20, 35),
         background = TRUE, map.polygon = NULL, raster.filename = raster_data,
         interpolation.model = FieldsKrigModel(10), col.palette = CreatePalette())

# 使用ggplot2绘制地理分布图
map.polygon <- getMap(resolution = "high")
pl <- ggtess3Q(q.matrix, coordinates, map.polygon = map.polygon)
pl + geom_path(data = map.polygon, aes(x = long, y = lat, group = group)) +
  coord_equal() +
  geom_point(data = as.data.frame(coordinates), aes(x = V1, y = V2), size = 0.2) +
  xlab("Longitude") + ylab("Latitude") + theme_bw()

# 进行基因组选择性扫描
p.values <- pvalue(tess3.obj, K = 5)
# 绘制p值的直方图
hist(p.values, col = "lightblue")

# 使用Benjamini-Hochberg方法进行FDR校正，并绘制“曼哈顿”图
L <- length(p.values)
fdr.level <- 1e-4
w <- which(sort(p.values) < fdr.level * (1:L)/L)
candidates <- order(p.values)[w]
plot(p.values, main = "Manhattan plot", xlab = "Locus id", ylab = "-log10(P-values)",
     cex = .3, col = "grey")
points(candidates, -log10(p.values)[candidates], pch = 19, cex = .2, col = "blue")
