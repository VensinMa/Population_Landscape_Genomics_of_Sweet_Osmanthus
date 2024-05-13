#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for simulating Gradient Forest and three genetic offset metrics 
#Install the dependencies package
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
rm(list=ls())
library(gradientForest)#version 0.1-24
#Load pop allele frequency
pika_snp <- read_table("pika.alfreq", col_names = TRUE)
#Load bio-climatic variables, in this study, we used bio5, bio6, and bio15
pika_env <- read.csv("pika_env.csv") %>% distinct()
#GF analysis
nSites <- dim(pika_snp)[1]
nSpecs <- dim(pika_snp)[2]
lev <- log2(nSites*0.368/2)
lev
gf <- gradientForest(cbind(pika_env,pika_snp), predictor.vars=colnames(pika_env), 
                     response.vars=colnames(pika_snp), ntree = 500,maxLevel = lev, 
                     corr.threshold = 0.5, compact = T, nbin = 201)
gf
#Predict the future, local genetic offset
#Load Qinghai-Tibet Plateau boundary, http://data.tpdc.ac.cn
QTP <- rgdal::readOGR("QTP.shp")
QTP_50k <- spsample(QTP,n=50000,"regular")#sample ~50000 grids
QTP_50k <- QTP_50k@coords %>% as.data.frame()
#Load current env layers
current_env <- raster("X")
QTP_50k_env <- raster::extract(current_env,QTP_50k,df=T)
vec <- c("bio5","bio6","bio15")
Trns_grid <- cbind(QTP_50k[1:2], 
                   predict(gf,QTP_50k_env[,vec]))
for (year in c(2070,2100)){
  for (i in c(126,585)){
    #future bio, CHELSA dataset: https://chelsa-climate.org/
    fut_bio <- list.files(paste0("D:\\Work\\High_ele_species\\Tibet_P_range/",
                                 year,"/SSP",i),
                          pattern = ".bil$", 
                          full.names = TRUE) %>% stack
    fut_bio_df <- raster::extract(fut_bio_both, 
                                  QTP_50k) %>% as.data.frame()
    Trns_grid_fut <- cbind(QTP_50k,predict(gf,fut_bio_df))
    Trns_grid_fut <- na.aggregate(Trns_grid_fut)#remove NA
    assign(paste0("Trns_grid_",year,"_SSP",i,sep = ""),Trns_grid_fut)
    df_go <- data.frame()
    for (j in 1:nrow(QTP_50k)){
     go_fut <- c(gv = rdist(Trns_grid[j,vec],Trns_grid_fut[j,vec]) )
    df_go <- bind_rows(df_go,go_fut)
    }
    assign(paste0("df_go_",year,"_SSP",i,sep = ""),df_gv)
  }
}
#Forward and reverse genetic offset, in this study, we only set suitable niche as target grids
require(parallel)
require(doParallel)
#Forward metric
predNames <- vec
popDatGF <- Trns_grid
popDatGF <- split(popDatGF, seq(nrow(popDatGF)))
#x:Future transformed genetic composition
#y:degrees of target grids
get_forwardoffset_suitable_each <- function(x,y){
  forwardOffsetGF <- foreach(i = 1:length(GTP_50K), .packages=c("fields","gdistance","geosphere")) %dopar%{
    
    #get the focal population
    onePopGF <- popDatGF[[i]]
    #get destination populations and add gf distance
    combinedDatGF <- y
    temp_trns <- left_join(y,x,by = c("Long", "Lat"))
    combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,predNames], temp_trns[,predNames]))
    ##Get metrics for the focal population
    #coordinate of focal population
    coordGF <- onePopGF[,c("Long","Lat")]
    #choose the pixels with the minimum fst
    minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
    #calculate the distance to the sites with minimum fst
    minCoordsGF["geoDist"] <- geosphere::distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
    #choose the pixels with the minimum geo
    minCoordsGF <- minCoordsGF[which(minCoordsGF$geoDist == min(minCoordsGF$geoDist)),]
    #if multiple sites have the same fst, and same least path cost, one is randomly chosen
    minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
    minPtGF <- minCoordsGF[,c("Long","Lat")]
    #write out 
    outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],  
               forwardOffset=minCoordsGF$gfOffset, 
               predDist=minCoordsGF$geoDist, 
               x2=minPtGF[[1]],y2=minPtGF[[2]])
    
  }
  forwardOffsetGF <- do.call(rbind, forwardOffsetGF)
}
#Backward metric
#x:future transformed genetic composition
#y:degrees of target grids
get_backwardoffset_suitable_each <- function(x,y){
  reverseOffsetGF <- foreach(i = 1:nrow(GTP_50K), .packages=c("fields","gdm","geosphere")) %dopar%{
    
    #get the focal population in future climate
    onePopGF <- x[i,]
    #get the focal population RGO based on suitable area
    temp <- y[,c("Long","Lat")]
    cur_stable_env <-left_join(temp,Trns_grid,by = c("Long", "Lat"))
    temp["gfOffset"] <- c(rdist(onePopGF[,predNames], cur_stable_env[,predNames]))
    #select the minimum one and calculate the distance
    temp_min <- temp[which(temp$gfOffset == min(temp$gfOffset)),]
    temp_min["dists"] <- geosphere::distGeo(p1=onePopGF[,1:2], p2=temp_min[,1:2])
    temp_min <- temp_min[which(temp_min$dists == min(temp_min$dists)),]
    temp_min <- temp_min[sample(1:nrow(temp_min),1),]
    #write out
    outGF <- c(x1=onePopGF[[1]], y1=onePopGF[[2]], reverseOffset=temp_min$gfOffset,
               predDist=temp_min$dists, x2=temp_min[[1]],y2=temp_min[[2]])
    
  }
  reverseOffsetGF <- do.call(rbind, reverseOffsetGF)
}