#!/usr/bin/env Rscript
library(xts)
library(bcp)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/BCP/TS.csv',index=0)

# Convert to numerics
#nums = as.numeric(myArgs)
#ts = ts(nums)

ser = as.numeric(ts.xts)
mvalue = bcp(ser, return.mcmc=TRUE)
#cat(cpts(mvalue))


write.csv(mvalue$posterior.prob,'./lib/BCP/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
