#!/usr/bin/env Rscript
library(xts)
library(cpm)

# max.R
# Fetch command line arguments
#myArgs <- commandArgs(trailingOnly = TRUE)
ts.xts <- read.zoo('./lib/CPM/TS.csv',index=0)


ser = as.numeric(ts.xts)

#vectors to hold the result
detectiontimes <- numeric()
changepoints <- numeric()
#use a Lepage CPM
cpm <- makeChangePointModel(cpmType="Lepage", ARL0=500)
i <- 0
while (i < length(ser)) {
      i <- i + 1
      #process each observation in turn
      cpm <- processObservation(cpm,ser[i])

      #if a change has been found, log it, and reset the CPM

      if (changeDetected(cpm) == TRUE) {
        #print(sprintf("Change detected at observation %d", i))
        detectiontimes <- c(detectiontimes,i)
        #the change point estimate is the maximum D_kt statistic
        Ds <- getStatistics(cpm)
        tau <- which.max(Ds)
        if (length(changepoints) > 0) {
          tau <- tau + changepoints[length(changepoints)]
        }
        changepoints <- c(changepoints,tau)
        #reset the CPM
        cpm <- cpmReset(cpm)
        #resume monitoring from the observation following the
        #change point
        i <- tau
} }
write.csv(changepoints,'./lib/CPM/CP.csv',row.names=FALSE, quote=FALSE, na="")
# cat will write the result to the stdout stream
#return cpts(mvalue)
