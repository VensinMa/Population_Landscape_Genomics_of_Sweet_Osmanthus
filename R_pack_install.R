.libPaths()

options(repos="http://mirrors.tuna.tsinghua.edu.cn/CRAN/")
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
install.packages("BiocManager")
packages=c("BiocManager", "dplyr", "DESeq2", "tidyverse", "ggplot2", "data.table", 
           "devtools", "shiny", "shinydashboard", "nadiv", "sommer", "agridat", 
           "agricolae", "desplot", "ggsci", "cowplot", "ggrepel", "ggsignif",
           "ggExtra", "ggstatsplot", "patchwork", "ggridges", "ggbeeswarm",
           "ggthemes", "gghalves", "ggvenn","ggcor", "ComplexHeatmap", "ggrepel", 
           "RIdeogram", "gggenes", "ggmsa",  "ggThemeAssist", "ggthemer", 
           "paletter", "PCAtools", "ggcor", "agricolae", "export","ggstasplot", 
           "customLayout", "writexl", "xlsx", "corrplot", "paletteer", "adegenet",
           "LEA","dismo", "gdistance", "factoextra", "FactoMineR", "dartR", "poppr",
           "SNPRelate", "treeio", "biomod2", "dismo", "ENMTools", "caret", "rinat",
           "bold", "taxize", "blockCV", "ggspatial", "tidyterra", "ggtext", "tidyterra",
           "maptools", "gdm", "usdm", "vegan", "ggthemr", "randomForest", "raster", 
           "caret", "geodata", "lfmm", "topGO", "Rgraphviz", "lavaan")
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    BiocManager::install(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
ipak(packages)
devtools::install_github('royfrancis/pophelper')
devtools::install_github(repo="nmatzke/BioGeoBEARS")
devtools::install_github("danlwarren/ENMTools") 

devtools::install_github("fitzLab-AL/gdm")
install.packages("conformal", repos="http://R-Forge.R-project.org")
install.packages("extendedForest", repos="http://R-Forge.R-project.org")
install.packages("gradientForest", repos="http://R-Forge.R-project.org")
