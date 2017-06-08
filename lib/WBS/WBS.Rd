#!/usr/bin/env Rscript
library(xts)
library(wbs)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/WBS/TS.csv',index=0)

# Convert to numerics
#nums = as.numeric(myArgs)
#ts = ts(nums)

ser = ts(ts.xts)
w <- wbs(ser)
w.cpt <- changepoints(w,penalty="mbic.penalty")
#cat(w.cpt)


write.csv(w.cpt$cpt.ic,'./lib/WBS/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
