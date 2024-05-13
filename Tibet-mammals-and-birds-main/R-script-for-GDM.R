#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for simulating Generalised Dissimilarity Modelling and three genetic offset metrics 
#Install the dependencies package
#install.packages("gdm")
rm(list = ls())
library(gdm)#version 1.4.2.2
#Load data
vec <- c("site", "Long","Lat","bio5","bio6","bio15")
env_tab <- read.csv("./pika_envTab.csv")[,vec]
pika_dis <- as.matrix(read.csv("pika_fst.csv",header = F))#Fst matrix
min <- min(pika_dis[upper.tri(pika_dis)])
max <- max(pika_dis)
pika_dis <- (pika_dis - min)/(max - min)
diag(pika_dis) <- 0
site <- env_tab$site %>% unique() 
pika_gdmdis  <- cbind(site,as.data.frame(pika_dis))
pika.gdmTab.dis <- formatsitepair(pika_gdmdis, bioFormat=3, XColumn="Long", YColumn="Lat",
                                  predData=env_tab, siteColumn="site")
pika_gdm <- gdm(pika.gdmTab.dis, geo=F)
#Predict the future, local genetic offset
#Load Qinghai-Tibet Plateau boundary, http://data.tpdc.ac.cn
QTP <- rgdal::readOGR("QTP.shp")
QTP_50k <- spsample(QTP,n=50000,"regular")#sample ~50000 grids
QTP_50k <- QTP_50k@coords %>% as.data.frame()
#Load current env layers
current_env <- raster("X")
QTP_50k_env <- raster::extract(current_env,QTP_50k,df=T)
for (year in c(2070,2100)){
  for (i in c(126,585)){
    #future bio, CHELSA dataset: https://chelsa-climate.org/
    fut_bio <- list.files(paste0("D:\\Work\\High_ele_species\\Tibet_P_range/",
                                 year,"/SSP",i),
                          pattern = ".bil$", 
                          full.names = TRUE) %>% stack
    fut_bio_both <- fut_bio_both[[c("bio5","bio6","bio15")]]%>%crop(buffer(QTP,1.5))%>%mask(TP)%>%stack
    fut_bio_df <- raster::extract(fut_bio_both, 
                                  QTP_50k) %>% as.data.frame()
    gdm_grid_fut <- cbind(QTP_50k,fut_bio_df)
    assign(paste0("gdm_grid_",year,"_SSP",i,sep = ""),gdm_grid_fut)
    fut.pred <- predict(pika_gdm, 
                       current_env, 
                            time = T,
                           predRasts=fut_bio_both)
    plot(fut.pred)
    assign(paste0("fut.pred_",year,"_SSP",i,sep = ""),fut.pred)
    fut.pred.df <- raster::extract(fut.pred,TP_50k)%>%as.data.frame
    assign(paste0("fut.pred.df_",year,"_SSP",i,sep = ""),fut.pred.df)
    
  }
}
#Forward and reverse genetic offset, in this study, we only set suitable niche as target grids
require(parallel)
require(doParallel)
#Forward metric
pred <- c("bio5","bio6","bio15")
popDat <- QTP_50k_env[c("Long","Lat","bio5","bio6","bio15")]
popDat <- data.frame(distance=1, weight=1, popDat)
popDat <- split(popDat, seq(nrow(popDat)))
removeIntercept <- function(mod,pred){
  adjust <- 0 - log(1-pred) - mod$intercept
  adjustDissim <- 1-exp(0-adjust)
  return(adjustDissim)
}
#x:degrees of target grids
#y:funture climatic conditions
get_forwardoffset_suitable_each <- function(x,y){
  #get target range
  temp_fut_df <- left_join(x,y)
  forwardOffsetGDM <- foreach(i = 1:nrow(QTP_50k), .packages=c("fields","gdistance","geosphere")) %dopar%{
    #get the focal population
    onePop <- popDat[[i]]
    #set up a dataframe where the first site is the focal population, and the second population
    setUp <- cbind(onePop,temp_fut_df,row.names = NULL)
    colnames(setUp) <- c("distance","weights",
                         "s1.xCoord", "s1.yCoord",paste("s1.", pred, sep=""), 
                         "s2.xCoord", "s2.yCoord",paste("s2.", pred, sep=""))
    
    #rearrange the colums for the gdm prediction
    dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                    paste("s1.", pred, sep=""), 
                    paste("s2.", pred, sep=""))]
    
    #do the prediction and set up a dataframe with second sites x/y and predicted Fst
    combinedDat <- predict(object=pika_gdm, dat, time=FALSE)
    combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(pika_gdm, combinedDat))
    
    ##Get metrics for the focal population
    #coordinate of focal population
    coord <- onePop[,c("Long","Lat")]
    #choose the pixels with the minimum fst and geo distance
    minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
    minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
    minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
    #if multiple sites have the same fst, and same distance, one is randomly chosen
    minCoords <- minCoords[sample(1:nrow(minCoords),1),]

    minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]

    #write out
    out <- c(x1=coord[[1]], y1=coord[[2]],
             forwardFst=minCoords$predFst, 
             predDist=minCoords$dists, 
             x2=minPt[[1]],y2=minPt[[2]])
  } 
  forwardOffsetGDM <- do.call(rbind, forwardOffsetGDM)
}
#
#Reverse
popDat_df <- QTP_50k_env[c("Long","Lat","bio5","bio6","bio15")]
#x:degrees of target grids
#y:funture climatic conditions
get_backwardoffset_suitable_each <- function(x,y){
  temp_curr_df <- left_join(x,popDat_df)
  reverseOffsetGDM <- foreach(i = 1:nrow(TP_50k), .packages=c("fields","gdm","geosphere")) %dopar%{
    #get the future focal population
    onePop <- y[i,]
    #set up a dataframe where the first site is the focal population, and the second population
    #are sites across the range
    setUp <- cbind(onePop,temp_curr_df,row.names = NULL)
    setUp <- data.frame(distance=1, weight=1, setUp)
    colnames(setUp) <- c("distance","weights",
                         "s1.xCoord", "s1.yCoord",paste("s1.", pred, sep=""), 
                         "s2.xCoord", "s2.yCoord",paste("s2.", pred, sep=""))
    #rearrange the colums for the gdm prediction
    dat <- setUp[,c("distance","weights","s1.xCoord", "s1.yCoord","s2.xCoord", "s2.yCoord",
                    paste("s1.", pred, sep=""), 
                    paste("s2.", pred, sep=""))]
    #do the prediction and set up a dataframe with second sites x/y and predicted Fst
    combinedDat <- predict(object=pika_gdm, dat, time=FALSE)
    combinedDat <- data.frame(dat[,c("s2.xCoord","s2.yCoord")], predFst=removeIntercept(pika_gdm, combinedDat))
    ##Get metrics for the focal population
    #coordinate of focal population
    coord <- onePop[,c("Long","Lat")]
    #choose the pixels with the minimum fst
    minCoords <- combinedDat[which(combinedDat$predFst == min(combinedDat$predFst)),]
    #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
    minCoords["dists"] <- distGeo(p1=coord, p2=minCoords[,1:2])
    minCoords <- minCoords[which(minCoords$dists == min(minCoords$dists)),]
    #if multiple sites have the same fst, and same distance, one is randomly chosen
    minCoords <- minCoords[sample(1:nrow(minCoords),1),]
    minPt <- minCoords[,c("s2.xCoord", "s2.yCoord")]
    #write out
    out <- c(x1=coord[[1]], y1=coord[[2]],
             reverseFst=minCoords$predFst, predDist=minCoords$dists, 
             x2=minPt[[1]],y2=minPt[[2]])
    
  }
  reverseOffsetGDM <- do.call(rbind, reverseOffsetGDM)
}