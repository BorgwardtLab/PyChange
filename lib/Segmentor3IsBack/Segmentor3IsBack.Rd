#!/usr/bin/env Rscript
library(xts)
library(Segmentor3IsBack)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/Segmentor3IsBack/TS.csv',index=0)


ser = as.numeric(ts.xts)
res = Segmentor(data=ser,model=3,Kmax=20,keep=TRUE)
Cr<-SelectModel(res,penalty='oracle',keep=FALSE)
Best<-BestSegmentation(res,K=Cr)


write.csv(getBreaks(res)[Cr,1:(Cr-1)],'./lib/Segmentor3IsBack/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
