#!/usr/bin/env Rscript
library(xts)
library(stepR)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/SMUCE/TS.csv',index=0)

# Convert to numerics
#nums = as.numeric(myArgs)
#ts = ts(nums)

ser = ts(ts.xts)
fit <- smuceR(ser, 1:length(ser), family="gauss")
#cat(w.cpt)


write.csv(fit$leftEnd,'./lib/SMUCE/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
