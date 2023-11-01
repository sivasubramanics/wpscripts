#!/usr/bin/env Rscript

args<-commandArgs(TRUE)
# Loading the Library
# if(!require(gplots)) install.packages("gplots")
# if(!require(LDheatmap)) install.packages("LDheatmap")
# if(!require(genetics)) install.packages("genetics")
# if(!require(ape)) install.packages("ape")
# if(!require(compiler)) install.packages("compiler")
# if(!require(grid)) install.packages("grid")
# if(!require(bigmemory)) install.packages("bigmemory")
# if(!require(EMMREML)) install.packages("EMMREML")
# if(!require(scatterplot3d)) install.packages("scatterplot3d")
# if(!require(lme4)) install.packages("lme4")


# library(multtest)
# library(gplots)
# library(LDheatmap)
# library(genetics)
# library(ape)
# library(EMMREML)
# library(compiler) #this library is already installed in R
# library("scatterplot3d")
# source("http://zzlab.net/GAPIT/gapit_functions.txt")
source("/lustre/BIF/nobackup/selva001/work/wp5/deepvariant/lsat/bremia_gwas/gapit_functions.txt")
# source("/data/Programs/GAPIT3/R/GAPIT.library.R")
# Loading the Gapit functions
# source("/home/sivasubramani/Programs/GAPIT/gapit_functions.txt")
# Read Trait data
# myY <- read.table(file="GWAS_Traits.txt", head = TRUE)
# Read HapMap file
# myG <- read.table(file="GWAS_Genotype.hmp.txt", head = FALSE)
# source("/data/Programs/GAPIT3/R/GAPIT.Multiple.Manhattan.R")
# METHOD <- args[1]
# source("http://zzlab.net/GAPIT/emma.txt")
# source("/data/Chickpea_136/GWAS/emma.txt")
# source("/home/sivasubramani/Programs/GAPIT/emma.txt")
# library(GAPIT)

CWD <- args[1]
PHE <- args[2]
HMP <- args[3]
inPCA <- 0 + as.numeric(args[4])
model <- args[5]
myY  <- read.table(PHE, head = TRUE)
myG <- read.delim(HMP, head = FALSE, sep="\t")

dir.create(CWD)
setwd(CWD)


cat ("PCA:", inPCA)

myGAPIT <- GAPIT(
  Y=myY,
  G=myG,
#   model="GLM",
  # model=c("GLM","MLM","Blink","MLMM","FarmCPU","SUPER"),# choose model
  model=model,
  # model=c("GLM"),# choose model
  PCA.total=inPCA,                                          # set total PCAs
#   cutOff=0.8,
  Inter.Plot=FALSE,                                      # perform interactive plot
  Multiple_analysis=TRUE,                               # perform multiple analysis
#   PCA.3d=TRUE,                                          # plot 3d interactive PCA
  file.output=TRUE
)




# # loading packages for GAPIT and GAPIT functions
# source("http://www.zzlab.net/GAPIT/GAPIT.library.R")
# source("http://www.zzlab.net/GAPIT/gapit_functions.txt")
#   # loading data set
# myY=read.table(file="http://zzlab.net/GAPIT/data/mdp_traits.txt", head = TRUE)
# myGD=read.table("http://zzlab.net/GAPIT/data/mdp_numeric.txt",head=T)
# myGM=read.table("http://zzlab.net/GAPIT/data/mdp_SNP_information.txt",head=T)
#   #myG=read.table(file="http://zzlab.net/GAPIT/data/mdp_genotype_test.hmp.txt", head = FALSE)
#   # performing simulation phenotype
# set.seed(198521)
# Para=list(h2=0.7,NQTN=20)
# mysimulation<-GAPIT(Para=Para,GD=myGD,GM=myGM)
# myY=mysimulation$Y


# myGAPIT <- GAPIT(
#   Y=myY[,c(1,2)],
#   GD=myGD,
#   GM=myGM,
#   model=c("GLM","MLM","SUPER","MLMM","FarmCPU","Blink"),# choose model
#   #model=c("FarmCPU"),
#   PCA.total=3,                                          # set total PCAs
#   NJtree.group=4,                                       # set the number of clusting group in Njtree plot
#   QTN.position=mysimulation$QTN.position,
#   Inter.Plot=TRUE,                                      # perform interactive plot
#   Multiple_analysis=TRUE,                               # perform multiple analysis
#   PCA.3d=TRUE,                                          # plot 3d interactive PCA
#   file.output=T
# )
