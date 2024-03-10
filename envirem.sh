# 加载必要的R包
library(terra)       # 用于光栅数据操作
library(sf)          # 空间矢量数据处理
library(rnaturalearth) # 获取地理数据
library(envirem)     # 环境变量计算
library(viridis)     # 高质量的颜色方案

# 定义西美地区裁剪区域的多边形
polyWKT <- "POLYGON((
  -130 50,
  -100 50,
  -100 30,
  -130 30,
  -130 50
))"
# 将WKT字符串转换为sf多边形对象
poly <- st_as_sfc(polyWKT, crs = 4326)

# 读取并准备输入光栅数据
inputDir <- "./chelsaClim/global"
files <- list.files(inputDir, pattern = '.tif$', full.names = TRUE)
clim <- rast(files) # 将所有光栅文件读入为一个堆栈

# 裁剪并掩膜光栅数据
clim <- crop(clim, poly) # 裁剪到多边形范围
land <- ne_download(scale = 50, type = 'land', category = 'physical', returnclass = 'sf')
clim <- mask(clim, vect(land)) # 使用下载的陆地多边形掩膜

# 聚合光栅数据以降低分辨率（可选）
clim <- aggregate(clim, factor = 4, fun = mean)

# 定义文件命名方案
assignNames(tmax = "CHELSA_tasmax_##_1981-2010_V.2.1",
            tmin = "CHELSA_tasmin_##_1981-2010_V.2.1",
            tmean = "CHELSA_tas_##_1981-2010_V.2.1",
            precip = "CHELSA_pr_##_1981-2010_V.2.1")

# 生成太阳辐射数据
# 注意：这一步假设1981-2010年间以1995年作为代表年份
ETsolradRasters(rasterTemplate = clim[[1]], year = 45, outputDir = "./chelsaClim/westernNA", overwrite = TRUE)
solrad <- rast(list.files('./chelsaClim/westernNA', pattern = 'solrad', full.names = TRUE))

# 计算ENVIREM变量
allGrids <- generateEnvirem(masterstack = clim, solradstack = solrad, var = 'all')

# 可视化并保存结果
outdir <- './chelsaClim/westernNA_envirem'
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# 遍历和可视化所有ENVIREM变量
par(mfrow = c(5, 4)) # 调整以适应您的变量数量和屏幕尺寸
for (i in 1:nlyr(allGrids)) {
  plot(allGrids[[i]], col = inferno(100), box = FALSE, axes = FALSE)
  title(main = names(allGrids)[i], cex = 0.65)
}

# 构建文件名并保存光栅数据到文件
filenames <- paste0(outdir, '/', names(allGrids), '.tif')
writeRaster(allGrids, filenames, overwrite = TRUE)
