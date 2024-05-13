#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for GEBVs.
#install.packages("rrBLUP")
library(rrBLUP) #v4.6.1
#Phenotype
pika_env <- read.csv("pika_env.csv")[c("bio5","bio6","bio15")]
pika_env$gid <- 1:nrow(pika_env)
#Genotype
pika.geno <- read_delim("./pika.genotype", 
                        "\t", col_names = FALSE, escape_double = FALSE)
#Calculate kinship matrix
pika.geno <- pika.geno - 1#convert 0/1/2 to -1/0/1
rownames(pika.geno) <- 1:nrow(pika_env)
A <- A.mat(pika.geno)
#Calculate GEBVs
#Using all three bio-climatic variables as input
pika_env_scale <- scale(pika_env)
y <- c(pika_env_scale[1:66,1],pika_env_scale[1:66,2],pika_env_scale[1:66,3])
env <- c(rep(1,66),rep(2,66),rep(3,66))
gid <- c(1:66,1:66,1:66)
data <- data.frame(y=y,env=env,gid=gid)
ans.all <- kin.blup(data,K=A,geno="gid",pheno="y",fixed="env")