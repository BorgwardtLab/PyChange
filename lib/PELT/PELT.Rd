#!/usr/bin/env Rscript
library(xts)
library(changepoint)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/PELT/TS.csv',index=0)

# Convert to numerics
#nums = as.numeric(myArgs)
#ts = ts(nums)

ser = ts(ts.xts)
mvalue = cpt.meanvar(ser, method='PELT')#,penalty='Manual',pen.value='0.05')
#cat(cpts(mvalue))


write.csv(cpts(mvalue),'./lib/PELT/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
