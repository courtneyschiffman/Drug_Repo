rm(list=ls())

# load results files
resultsFiles <- grep("p1", list.files("./data"), value=T)
for(x in resultsFiles) { load(paste0("./data/", x)) }

plotResults <- function(permResults) {
  nSim <- length(permResults)
  nStats <- length(permResults[[1]]$simStats)
  statsNames <- colnames(permResults[[1]]$permStats)
  par(mfcol=c(nStats, nSim))
  # compute xlim for plots
  xMax <- apply(sapply(permResults,
                       function(x) {
                         apply(x$permStats, 2, max)
                       }), 1, max)
  xMin <- apply(sapply(permResults,
                       function(x) {
                         apply(x$permStats, 2, min)
                       }), 1, min)
  xLims <- matrix(NA, nrow=nStats, ncol=2)
  simStats <- sapply(permResults, function(x) x$simStats)
  for(i in 1:nStats) {
    xLims[i,1] <- min(xMin[i], simStats[i,]) - 0.1
    xLims[i,2] <- max(xMax[i], simStats[i,]) + 0.1
  }
  for(i in 1:nSim) {
    for(j in 1:nStats) {
      hist(permResults[[i]]$permStats[,j],
           breaks=seq(from=xLims[j,1], to=xLims[j,2], length.out=20),
           main=paste0("simulation", i, ": ", statsNames[j]),
           xlab=paste("p-value:", permResults[[i]]$pvalues[j]))
      abline(v=permResults[[i]]$simStats[j], col="red")
    }
  }
  par(mfrow=c(1,1))
}

plotResults(p1_0.2_p2_0.05_p3_0)
plotResults(p1_0.2_p2_0.05_p3_0.5)
plotResults(p1_0.2_p2_0.05_p3_1)
