rm(list=ls())

# load results files
resultsFiles <- grep("noidr", list.files("./data"), value=T)
for(x in resultsFiles) { load(paste0("./data/", x)) }

plotResults <- function(permResults) {
  require(ggplot2)
  require(reshape)
  require(gridExtra)
  nSim <- length(permResults)
  nPerm <- nrow(permResults[[1]]$permStats)
  nStats <- ncol(permResults[[1]]$permStats)
  statsNames <- colnames(permResults[[1]]$permStats)
  # compute average p-value over simulations
  allPvalues <- as.data.frame(t(sapply(permResults, function(x) x$pvalues)))
  colnames(allPvalues) <- statsNames
  avgPvalues <- colMeans(allPvalues)
  # collect stats over all simulations
  allStats <- as.data.frame(t(sapply(permResults, function(x) x$simStats)))
  colnames(allStats) <- paste0("simulation", statsNames)
  # collect stats over all permutations
  allPerms <- data.frame(matrix(NA, nrow=nSim*nPerm, ncol=nStats))
  for(i in 1:nSim) {
    allPerms[((i-1)*1000+1):(nPerm*i),] <- permResults[[i]]$permStats
  }
  colnames(allPerms) <- paste0("permutation", statsNames)
  # combine simulation and permutation stats
  facetTitles <- paste0(statsNames, ", p=", avgPvalues)
  names(facetTitles) <- statsNames
  allDat <- rbind(melt(allPerms), melt(allStats))
  allDat$stat <- sub("permutation", "", sub("simulation", "", allDat$variable))
  allDat$stat <- facetTitles[allDat$stat]
  if(nStats == 3) {
    allDat$type <- sub("gsea", "", sub("spearman", "", sub("trunSpearman", "", allDat$variable)))
  }
  if(nStats ==4) {
    allDat$type <- sub("gsea", "", sub("spearman", "", sub("trunSpearman", "", sub("IDR.func", "", allDat$variable))))
  }
  # plot results
  params <- unlist(strsplit(deparse(substitute(permResults)), split="_"))
  # plotTitle <- paste(sapply(1:(length(params)/2),
  #                     function(x) {
  #                       paste(params[2*(x-1)+1], params[2*x], sep="=")
  #                     }),
  #                    collapse=", ")
  plotTitle <- paste(100*as.numeric(params[6]), "% same direction,",
                     100*as.numeric(params[8]), "% diff. direction")
  histPlot <- ggplot(subset(allDat, allDat$stat==x), aes(value)) +
    geom_histogram(data=subset(allDat, allDat$type=="permutation"),
                   binwidth=0.05, fill="blue", alpha=0.5) +
    geom_dotplot(data=subset(allDat, allDat$type=="simulation"),
                 binwidth=0.05, fill="red", alpha=0.5) +
    ggtitle(plotTitle) + xlab("") +
    facet_grid(.~stat) + theme_bw()
  if(nStats == 3) {
    pvaluePlot <- ggplot(melt(allPvalues), aes(variable, value)) +
      geom_violin(fill="green", alpha=0.2) + xlab("") + ylab("pvalues") +
      theme_bw()
  }
  if(nStats ==4) {
    pvaluePlot <- ggplot(melt(allPvalues), aes(variable, value)) +
    geom_dotplot(binaxis="y", stackdir="center") + xlab("") + ylab("pvalues") +
    theme_bw()
  }
  grid.arrange(histPlot, pvaluePlot, nrow=2)
}

plotResults(p1_0.05_p2_0.05_p3_0_p4_0)
plotResults(p1_0.05_p2_0.05_p3_0_p4_0.5)
plotResults(p1_0.05_p2_0.05_p3_0.5_p4_0)
plotResults(p1_0.05_p2_0.05_p3_0.5_p4_0.5)
plotResults(p1_0.05_p2_0.05_p3_0_p4_1)
plotResults(p1_0.05_p2_0.05_p3_1_p4_0)

## new results w/ idr
newResults <- grep("noidr", grep("p4", list.files("./data"), value=T), invert=T, value=T)
for(x in newResults) { load(paste0("./data/", x)) }
plotResults(p1_0.05_p2_0.05_p3_0_p4_0)
plotResults(p1_0.05_p2_0.05_p3_0.5_p4_0)
plotResults(p1_0.05_p2_0.05_p3_1_p4_0)
