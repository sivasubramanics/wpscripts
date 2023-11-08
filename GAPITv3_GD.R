#!/usr/bin/env Rscript

# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT3",force=TRUE)
library(GAPIT)

args<-commandArgs(TRUE)

CWD <- args[1] # Ouput directory
PHE <- args[2] # Phenotype file
GD <- args[3] # Genotype data
GM <- args[4] # Genotype map
MODEL <- args[5] # Model

# load genotype data and phenotype data
myGD=read.table(GD,header=T,sep="\t")
myGM=read.table(GM,header=T,sep="\t")
myY=read.table(PHE,header=T,sep="\t")

# create a directory for GAPIT output
dir.create(CWD)
setwd(CWD)

# run GAPIT
myGAPIT=GAPIT(
  Y=myY, #fist column is ID
  GD=myGD,
  GM=myGM,
  PCA.total=3,
  model=MODEL,
  Multiple_analysis=TRUE)
