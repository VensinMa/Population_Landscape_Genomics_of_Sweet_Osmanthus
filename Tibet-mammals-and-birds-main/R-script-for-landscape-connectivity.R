#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for Landscape connectivity analysis.
#The first part shows how to use ResistanceGA optimizing resistance surfaces, and the second part used for MLPE
############
#first part#
############
#install.packages("devtools")
#devtools::install_github("wpeterman/ResistanceGA") # v4.2-2
library(ResistanceGA)
library(raster)
library(parallel)
library(doParallel)
#Before optimizing, you need to stack all environment layers into one raster file.
#For more details, plz refer to https://github.com/wpeterman/ResistanceGA/
#Load raster file
load("pika_resis_predata.Rdata")
#Optimization with CIRCUITSCAPE (v4) in julia platform
JULIA_HOME <- "/home/prunella/HIC_data_flow/phylonet/julia-1.7.0/julia-1.7.0/bin"
JuliaCall::julia_setup(JULIA_HOME)
#set threads number you want to use
cl <- makePSOCKcluster(5)
registerDoParallel(cl)
jl.inputs <- jl.prep(n.Pops = 12,
                     response = pika_fst[lower.tri(pika_fst)],
                     CS_Point.File = sp.dat,
                     JULIA_HOME = JULIA_HOME,
                     run_test = F)
GA.inputs <- GA.prep(ASCII.dir = "./resistance",
                     Results.dir = "./resistance/",
                     min.cat = 1,
                     max.cat = 100,
                     max.cont = 100,
                     #set surface transformations for time saving if you wish
                     #select.trans = list(c(2,8),c(2,8),c(2,8),c(1,3),c(5,7),c(2,8),c(1,3),NA),
                     method = "AIC",
                     seed = 555,
                     parallel=cl,
                     quiet = TRUE)
# Export info to cluster
clusterExport(cl=cl,varlist=c("jl.inputs","GA.inputs","sp.dat","pika_fst")) # list everything you call in ro GA.inputs and gdist
clusterEvalQ(cl=cl, .libPaths("/home/sparrow/.conda/envs/rgdal/lib/R/library")) # set path to where your R library is
clusterCall(cl=cl, library, package = "ResistanceGA", character.only = TRUE)

jl.optim <- SS_optim(jl.inputs = jl.inputs,
                     GA.inputs = GA.inputs)

stopCluster(cl)
#############
#second part#
#############
#install.packages("lme4")
#install.packages("MuMIn")
library(lme4) #v1.1-27.1
library(MuMIn) #1.43.1
#read fst matrix
pika_fst <- read.csv("pika_fst.csv",header = F)
#Geo distance resistance surface
geo_res <- read.csv("./resistance/Results/Distance_jlResistMat.csv",header = F)
geo_res <- as.data.frame(geo_res[ lower.tri(geo_res)])
#Climate resistance surface
bio5_res <- read.csv("./resistance/Results/bio5_resistance_jlResistMat.csv",header = F)
bio5_res <- bio5_res[ lower.tri(bio5_res)]%>% as.data.frame()

bio6_res <- read.csv("./resistance/Results/bio6_resistance_jlResistMat.csv",header = F)
bio6_res <- bio6_res[ lower.tri(bio6_res)]%>% as.data.frame()

bio15_res <- read.csv("./resistance/Results/bio15_resistance_jlResistMat.csv",header = F)
bio15_res <- bio15_res[ lower.tri(bio15_res)]%>% as.data.frame()
#Land use resistance surface
grass_res <- read.csv("./resistance/Results/grass_resistance_jlResistMat.csv",header = F)
grass_res <- grass_res[ lower.tri(grass_res)]%>% as.data.frame()

forest_res <- read.csv("./resistance/Results/forest_resistance_jlResistMat.csv",header = F)
forest_res <- forest_res[ lower.tri(forest_res)]%>% as.data.frame()
#Topology resistance surface
alt_res <- read.csv("./resistance/Results/alt_resistance_jlResistMat.csv",header = F)
alt_res <- alt_res[ lower.tri(alt_res)]%>% as.data.frame()

slope_res <- read.csv("./resistance/Results/slope_resistance_jlResistMat.csv",header = F)
slope_res <- slope_res[ lower.tri(slope_res)]%>% as.data.frame()

river_res <- read.csv("./resistance/Results/river_resistance_jlResistMat.csv",header = F)
river_res <- river_res[ lower.tri(river_res)]%>% as.data.frame()
#MLPE
dim <- read.csv("pika_MLPE_df.csv",header = T)#population effect matrix
df <- cbind(fst_dis,geo_res,bio5_res,bio6_res,bio15_res,
            grass_res,forest_res,
            alt_res,slope_res,river_res)
names(df) <- c("Genetic_distance","IBD","MTWM","MTCM","PS",
               "Grass","Forest",
               "Elevation","Slope","River")
names(df) <- c("Genetic_distance","IBD","MTWM","MTCM","PS",
               "Grass","Forest",
               "Elevation","Slope","River")
df <- cbind(dim,df)
str(df)
df <- df %>% mutate(across(c(3:12),scale))#standardization
#Define a function
MLPE <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = TRUE)
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}
#For comparing AIC
MLPEnoREML <- function(variables, data) {
  mod2 <- lme4::lFormula(variables, data = data, REML = FALSE)#REML = FALSE
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  mod_2 <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  mod2$reTrms$Zt <- ZZ
  
  # Refit the model
  dfun <- do.call(lme4::mkLmerDevfun, mod2)
  opt <- lme4::optimizeLmer(dfun)
  modelout <- lme4::mkMerMod(environment(dfun), opt, mod2$reTrms,fr = mod2$fr)
  return(modelout)
}
#model for compare
(model_lc_1_aic <- MLPEnoREML(Genetic_distance ~ IBD + (1|pop1),df))#geo hypothesis
(model_lc_2_aic <- MLPEnoREML(Genetic_distance ~ MTWM + MTCM + PS + (1|pop1),df))#climate hypothesis
(model_lc_3_aic <- MLPEnoREML(Genetic_distance ~ Grass + Forest + (1|pop1),df))#land use hypothesis
(model_lc_4_aic <- MLPEnoREML(Genetic_distance ~ MTWM + MTCM + PS + Grass + Forest + (1|pop1),df))#env hypothesis
(model_lc_5_aic <- MLPEnoREML(Genetic_distance ~ Elevation + River + Slope  + (1|pop1),df))#topology hypothesis
Models <- list(
                   IBD=model_lc_1_aic, 
                   GCC=model_lc_2_aic, 
                   LU=model_lc_3_aic,
                   GCC_LU=model_lc_4_aic,
                   Topo = model_lc_5_aic)
#Compute the Marginal R2 (the variance explained by fixed factors) and the Conditional R2 (variance explained by both fixed and random factors, whole model). 
df <- data.frame()
for (i in 1:5){
  df_modle <- MuMIn::r.squaredGLMM(Models[[i]]) %>% as.data.frame()
  df <- rbind(df,df_modle)
}
write.csv(df, "pika_resis_MLPE.csv")