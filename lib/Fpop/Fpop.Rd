#!/usr/bin/env Rscript
library(xts)
library(fpop)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/Fpop/TS.csv',index=0)


ser = as.numeric(ts.xts)
cp = Fpop(ser,length(ser)/20)

#print(cp)

write.csv(cp$t.est,'./lib/Fpop/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
