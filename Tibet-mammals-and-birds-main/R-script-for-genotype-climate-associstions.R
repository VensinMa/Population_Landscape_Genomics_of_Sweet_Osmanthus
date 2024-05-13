#Notes:as we have three species used in our analysis, for simplicity, this script will just use Ochotona curzoniae (pika) as an example.
#This script used for genotype-climate associations, including LFMM & RDA analysis.
############
####LFMM####
############
#install.packages("devtools")
#devtools::install_github("bcm-uga/LEA")
library(LEA)#version 3.4.0
#read genotype data
pika.geno <- read_delim("./pika.genotype", 
	"\t", col_names = FALSE, escape_double = FALSE) %>% as.matrix()
#read climate data
pika <- read.csv("pika_bio.csv",header = T) %>% as.matrix()
#run lfmm
bio <-c(5,6,15)
p_num<-c(0.05)
for(i in bio){
  for(p in p_num){
    
    mod<- lfmm2(input = pika.geno, env = pika[,i], K = 1)
    pv<- lfmm2.test(object = mod,
                    input = pika.geno,
                    env = pika[,i],
                    linear = TRUE)
    qval<- qvalue(pv$pvalues)$qvalues
    alpha <- p
    outliers<- which(qval< alpha)
    write.table(outliers,paste("pika_bio",i,"_k1_",p,".txt",sep=""),row.names = FALSE,col.names = FALSE)
  }
}
############
####RDA#####
############
#For more Redundancy Analysis (RDA) details, please refer to https://popgen.nescent.org/2018-03-27_RDA_GEA.html
library(vegan)    # Used to run RDA,v2.5-7
#run rda
pika.rda <- rda(pika.geno ~ ., data=pika[,c(5,6,15)], scale=T)
#Check for significance
signif.full <- anova.cca(pika.rda, parallel=getOption("mc.cores"))
signif.axis <- anova.cca(pika.rda, by="axis", parallel=getOption("mc.cores"))
#Identify candidate SNPs involved in local adaptation
load.rda <- scores(pika.rda, choices=c(1:3), display="species")#The choosen number depends on the numbers of significance axis
SNP_df <- as.data.frame(row.names(load.rda))
SNP_df$pos <- 1:nrow(SNP_df)
names(SNP_df) <- c("snp", "pos")
#3 standard deviation cutoff
outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
cand1 <- outliers(load.rda[,1],3) 
cand2 <- outliers(load.rda[,2],3) 
cand3 <- outliers(load.rda[,3],3)
ncand <- length(cand1) + length(cand2) + length(cand3)
cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))
cand2 <- cbind.data.frame(rep(2,times=length(cand2)), names(cand2), unname(cand2))
cand3 <- cbind.data.frame(rep(3,times=length(cand3)), names(cand3), unname(cand3))
colnames(cand1) <- colnames(cand2) <- colnames(cand3) <-c("axis","snp","loading")
cand <- rbind(cand1, cand2, cand3)
cand$snp <- as.character(cand$snp)
#correlations of each candidate SNP with the environmental predictors
foo <- matrix(nrow=nrow(cand), ncol=3)  # 3 columns for 3 predictors
colnames(foo) <- c("bio5","bio6","bio15")
for (i in 1:length(cand$snp)) {
  nam <- cand[i,2]
  snp.gen <- pika.geno[,nam]
  foo[i,] <- apply(pika[,c(5,6,15)],2,function(x) cor(x,snp.gen))
}
cand <- cbind.data.frame(cand,foo)
length(cand$snp[duplicated(cand$snp)])  
cand <- cand[!duplicated(cand$snp),] # remove duplicate detections
head(cand)
cand.pos <- left_join(cand, SNP_df, by = "snp")
#save outliers
write.table(cand.pos[,7], "./pika.cand.txt", quote = F, row.names = F,col.names = F)
