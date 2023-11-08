#!/usr/bin/env Rscript

# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT)

args<-commandArgs(TRUE)

CWD <- args[1]
PHE <- args[2]
GD <- args[3]
GM <- args[4]
MODEL <- args[5]

dir.create(CWD)
setwd(CWD)

myGD=read.table(GD,header=T,sep="\t")
myGM=read.table(GM,header=T,sep="\t")
myY=read.table(PHE,header=T,sep="\t")

myGAPIT=GAPIT(
  Y=myY, #fist column is ID
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  model=MODEL,
  Multiple_analysis=TRUE)
