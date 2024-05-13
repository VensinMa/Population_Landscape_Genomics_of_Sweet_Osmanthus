#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for Ecological niche modelling
#The first part shows how to use SDMtune tuning Maxent, and the second part used for ensemble modeling using biomod2
#Install the dependencies package
#install.packages("SDMtune")
#install.packages("devtools")
#devtools::install_github("HemingNM/ENMwizard")
#install.packages("biomod2", dependencies = TRUE)
rm(list = ls())
library(megaSDM)#version 1.0.1
library(SDMtune)#version 1.1.5
library(biomod2)#version 3.5.1
############
#first part#
############
#load current env data
QTP <- rgdal::readOGR("QTP.shp")
input_TA <- list.files("D://Work//High_ele_species//Tibet_P_range/current/",
                       pattern = ".bil$", 
                       full.names = TRUE) %>%stack
landu_TA <- list.files("D://Work//High_ele_species//Tibet_P_range/LUH2/",
                       pattern = "current", 
                       full.names = TRUE) %>%stack
names(landu_TA) <- c("forest","grass")
current_env <- stack(input_TA,landu_TA)%>%crop(buffer(QTP,1.5))%>%
  mask(buffer(QTP,1.5))%>%stack
#set output dir
envoutput <- "TestRun"
#set extent
ex <- extent(c(70, 110, 20, 45))
spplist <- c("Ochotona curzoniae")#pika locations
occ_output <- "occurrences"
#extract bio value
occlist <- list.files(occ_output, pattern = ".csv", full.names = TRUE)
OccurrenceManagement(occlist = occlist,
                     output = occ_output,
                     envextract = TRUE,
                     envsample = F,
                     nbins = 25,
                     envdata = current_env)
#backgroud point
buff_output <- "TestRun/buffers"
# Generates buffers for each species.
BackgroundBuffers(occlist = occlist,
                  envdata = current_env,
                  buff_output,
                  buff_distance = 1,
                  ncores = 2)
nbg <- 10000#backgroud
spatial_weights <- 0.5
sampleMethod <- "Varela"#subsampling method
bufflist <- list.files(buff_output, pattern = ".shp$", full.names = TRUE)
bg_output <- "TestRun/backgrounds"
BackgroundPoints(spplist = spplist,
                 envdata = current_env,
                 output = "TestRun/backgrounds",
                 nbg = nbg,
                 spatial_weights = spatial_weights,
                 buffers = bufflist,
                 method = sampleMethod,
                 ncores = 4)
#select Maxent hyperparameter value using SDMturn
pika_sdmturn <- prepareSWD(species = "Ochotona curzoniae", 
                           p = read.csv(occlist[1])[,2:3], 
                           a = read.csv(bglist[1])[,2:3], 
                           env = current_env[[c("bio5","bio6","bio15","grass","forest")]])
datasets <- trainValTest(pika_sdmturn, test = 0.3, only_presence = TRUE)
train <- datasets[[1]]
test <- datasets[[2]]
pika_sdmturn_t <- train("Maxent",data = train)
h <- list(reg = seq(0.5, 5, 0.5), 
          fc = c("l", "p", "q", "h", "lp", "lq", "lh", "pq", "ph", "qh", 
                 "lpq", "lph", "lqh","pqh", "lpqh"))
output <- gridSearch(pika_sdmturn_t, hypers = h, 
                     metric = "aicc", env = current_env[[c("bio5","bio6","bio15","grass","forest")]],
                     test = test,save_models = FALSE)
head(output@results[order(output@results$delta_AICc), ])#check the optimal one

#############
#second part#
#############
#load bio-interaction variable
finch_occur_surface <- raster("X.tif")
current_env <- stack(current_env,finch_occur_surface)
#set para for biomod2
DataSpecies <- read.csv(occlist[1])[,2:3]
bg_df <- read.csv(bglist[1])[,2:3]
pa_random <- rbind(DataSpecies,bg_df)
pa_random$PA1 <- as.logical("TRUE")
myRespName <- 'pika'
myResp <- as.numeric(c(rep(1,nrow(DataSpecies)),rep(NA,nrow(bg_df))))
#,"interact"
myExpl <- current_env[[c("bio5","bio6","bio15","grass","forest","interact")]]%>%stack#c("bio5","bio6","bio15","grass","forest")
myExpl_proj <-  current_env[[c("bio5","bio6","bio15","grass","forest","interact")]]%>%mask(TP)%>%stack

#data formating, selecte PA randomly
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl ,
                                     resp.xy =pa_random[,1:2],
                                     resp.name = myRespName,
                                     PA.strategy = 'user.defined',
                                     PA.table = pa_random[,3]%>%as.matrix()
)
#set maxent parameter based on the result of model tuning
myBiomodOption_2 <- BIOMOD_ModelingOptions(MAXENT.Phillips = list( path_to_maxent.jar = paste(system.file(package="dismo"), "/java", sep=''),
                                                                   memory_allocated = 4096,
                                                                   background_data_dir = 'default',
                                                                   maximumbackground = 'default',
                                                                   maximumiterations = 200,
                                                                   visible = FALSE,
                                                                   linear = F,
                                                                   quadratic = F,
                                                                   product = T,
                                                                   threshold = FALSE,
                                                                   hinge = T,
                                                                   lq2lqptthreshold = 80,
                                                                   l2lqthreshold = 10,
                                                                   hingethreshold = 15,
                                                                   beta_threshold = -1,
                                                                   beta_categorical = -1,
                                                                   beta_lqp = -1,
                                                                   beta_hinge = -1,
                                                                   betamultiplier = 3.0,
                                                                   defaultprevalence = 0.5))
#
set.seed(123)
myBiomodModelOut <- BIOMOD_Modeling( myBiomodData, 
                                     models = c('MAXENT.Phillips','GBM','RF','MARS'), 
                                     models.options = myBiomodOption_2, 
                                     NbRunEval=5, 
                                     DataSplit=70, 
                                     Yweights=NULL, 
                                     VarImport=3, 
                                     models.eval.meth = c('TSS','ROC'),
                                     SaveObj = TRUE,
                                     rescal.all.models = T,
                                     do.full.models = FALSE,
                                     modeling.id='pika')
model.eva<- get_evaluations(myBiomodModelOut)
model.eva["TSS","Testing.data",,,]
model.eva["ROC","Testing.data",,,]
models_scores_graph(myBiomodModelOut, by = "models" ,
                    metrics = c("TSS","ROC"))
var_imp <- get_variables_importance(myBiomodModelOut)
write.csv(var_imp,file=paste(myRespName,"variable_importance.csv",sep="_"),row.names=T)
#EnsembleModeling
myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                      chosen.models = 'all',
                                      em.by='all',
                                      eval.metric = c('TSS','ROC'),
                                      eval.metric.quality.threshold = c(0.6,0.8),
                                      models.eval.meth = c('TSS','ROC'),
                                      prob.mean = T,
                                      prob.cv = F,
                                      prob.ci = F,
                                      prob.ci.alpha = 0.05,
                                      prob.median = T,
                                      committee.averaging = F,
                                      prob.mean.weight = T,
                                      prob.mean.weight.decay = 'proportional' )
model.eva.EM <- get_evaluations(myBiomodEM)
model.eva.EM
#project
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl_proj,
  proj.name = 'current', 
  selected.models = 'all', 
  binary.meth = 'TSS', 
  compress = 'gzip', 
  clamping.mask = T, 
  output.format = '.grd',
  do.stack=T)
mod_projPres <- get_predictions(myBiomodProj)
#project EnsembleModeling
myBiomodEF <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                         projection.output = myBiomodProj,
                                         binary.meth ="TSS",
                                         proj.name = 'current_EM')
projPresEnsemble <- get_predictions(myBiomodEF)
#project to future

for (year in c(2070,2100)){
  for (i in c(126,585)){
    #future bio
    fut_bio <- list.files(paste0("D:\\Work\\High_ele_species\\Tibet_P_range/",
                                 year,"/SSP",i),
                          pattern = ".bil$", 
                          full.names = TRUE) %>% stack #%>%crop(TP) %>% mask(TP)
    fut_lu <- list.files("D://Work//High_ele_species//Tibet_P_range/LUH2/",
                         pattern = paste0(i,"_",year), 
                         full.names = TRUE) %>% stack #%>%crop(TP) %>% mask(TP)
    names(fut_lu) <- c("forest","grass")
    fut_bio_both <- stack(fut_bio,fut_lu)
    fut_bio_both <- fut_bio_both[[c("bio5","bio6","bio15","grass","forest")]]%>%crop(buffer(TP,1.5))%>%mask(TP)%>%stack
    fut_bio_both <- stack(fut_bio_both,myExpl_proj[["interact"]])
    
    #project
    myBiomodProj_loop <- BIOMOD_Projection(
      modeling.output = myBiomodModelOut,
      new.env = fut_bio_both,
      proj.name = paste0(year,"_SSP",i), 
      selected.models = 'all', 
      binary.meth = 'TSS', 
      compress = 'gzip', 
      clamping.mask = T, 
      output.format = '.grd',
      do.stack=T)
    myBiomodEF_loop <- BIOMOD_EnsembleForecasting(EM.output = myBiomodEM,
                                                  projection.output = myBiomodProj_loop,
                                                  binary.meth ="TSS",
                                                  proj.name = paste0(year,"_SSP",i,"_EM"))
  }
}
