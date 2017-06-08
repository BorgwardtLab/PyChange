#!/usr/bin/env Rscript
library(xts)
library(ecp)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/E_Divise/TS.csv',index=0)

# Convert to numerics
#nums = as.numeric(myArgs)
#ts = ts(nums)

ser = ts(ts.xts)
ser = matrix(c(ser), ncol = 1)
w.cpt <- e.divisive(ser)$estimates
#cat(w.cpt)


write.csv(w.cpt,'./lib/E_Divise/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
