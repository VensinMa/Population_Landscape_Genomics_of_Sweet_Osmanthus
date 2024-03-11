# 加载必要的R包
install.packages("envirem")
library(terra)       # 用于光栅数据操作
library(sf)          # 空间矢量数据处理
library(rnaturalearth) # 获取地理数据
library(envirem)     # 环境变量计算
library(viridis)     # 高质量的颜色方案

# 定义中国裁剪区域的多边形
polyWKT <- "POLYGON((96 33, 124 33, 124 18, 96 18, 96 33))"

# 将WKT字符串转换为sf多边形对象
poly <- st_as_sfc(polyWKT, crs = 4326)
poly

land <- ne_download(scale = 10, type = 'land', category = 'physical', returnclass = 'sf')
land <- st_geometry(land)
plot(poly, border = 'blue', xpd = NA)
plot(land, border = NA, col = gray(0.9, alpha = 0.5), add = TRUE, xpd = NA)

# 读取并准备输入光栅数据
inputDir <- "E:/mwx/envirem"
files <- list.files(inputDir, pattern = '.tif$', full.names = TRUE)
clim <- rast(files) # 将所有光栅文件读入为一个堆栈

verifyRasterNames(clim)
# 裁剪并掩膜光栅数据
clim <- crop(clim, poly) # 裁剪到多边形范围
land <- ne_download(scale = 50, type = 'land', category = 'physical', returnclass = 'sf')
clim <- mask(clim, vect(land)) # 使用下载的陆地多边形掩膜

# 聚合光栅数据以降低分辨率（可选）
clim <- aggregate(clim, factor = 4, fun = mean)
# ?aggregate
# 定义文件命名方案
assignNames(tmax = "CHELSA_tmax_##_2041-2070_SSP126",
            tmin = "CHELSA_tmin_##_2041-2070_SSP126",
            tmean = "CHELSA_tmean_##_2041-2070_SSP126",
            precip = "CHELSA_precip_##_2041-2070_SSP126")

# 生成太阳辐射数据
# 注意：这一步假设1981-2010年间以1995年作为代表年份 1950+
ETsolradRasters(rasterTemplate = clim[[1]], year = 105, outputDir = "E:/mwx/envirem/test", overwrite = TRUE)
solrad <- rast(list.files('E:/mwx/envirem/test', pattern = 'solrad', full.names = TRUE))

# 计算ENVIREM变量
allGrids <- generateEnvirem(masterstack = clim, solradstack = solrad, var = 'all')

# 可视化并保存结果
outdir <- 'E:/mwx/envirem/test_envirem'
if (!dir.exists(outdir)) {
  dir.create(outdir)
}

# 遍历和可视化所有ENVIREM变量
par(mfrow = c(5, 4)) # 调整以适应您的变量数量和屏幕尺寸
for (i in 1:nlyr(allGrids)) {
  plot(allGrids[[i]], col = inferno(100), box = FALSE, axes = FALSE)
  title(main = names(allGrids)[i], cex = 0.65, line = 3)
}

# 构建文件名并保存光栅数据到文件
filenames <- paste0(outdir, '/', names(allGrids), '.tif')
writeRaster(allGrids, filenames, overwrite = TRUE)



