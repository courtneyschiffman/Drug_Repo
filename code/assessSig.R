
assessSig <- function(expDat, simFuns, simFunsDir, nPerm, deg=49) {
  ###  computes similarity scores for nPerm permutations
  # expDat <- expression data: data.frame with x = signature1 (disease) and y=signature2 (drug)
  # simFuns: vector of similarity functions
  # simFunsDir: logical vector of whether more positive values mean better profile matching
  # nPerm: number of permutations
  simStats <- sapply(simFuns,
                     function(fun) {
                       do.call(fun, args=list(x=expDat$x,
                                              y=expDat$y,
                                              deg=deg))
                     })
  permStats <- foreach(iterators::icount(nPerm), .combine='rbind') %dopar%
    sapply(simFuns,
           function(fun) {
             do.call(fun, args=list(x=sample(expDat$x),
                                    y=expDat$y,
                                    deg=deg))
           })
  pvalues <- sapply(1:length(simFuns),
                    function(i) ifelse(simFunsDir[i],
                                       mean(permStats[,i]>simStats[i]),
                                       mean(permStats[,i]<simStats[i])))
  return(list(simStats=simStats,
              permStats=permStats,
              pvalues=pvalues))
}

## example
# simDat <- data.frame(x=rt(20000, df=49)+log(5), y=rt(20000, df=49))
# spearman <- function(x, y, deg) {
#   cor(x, y, method="spearman")
# }
# trunSpearman <- function(x, y, deg) {
#   x.pvalues <- pt(x, df=deg)
#   x.pvalues[x.pvalues>0.5] <- 1-x.pvalues[x.pvalues>0.5]
#   x.pvalues <- 2*x.pvalues
#   use <- p.adjust(x.pvalues, method="fdr")<0.05
#   cor(x[use], y[use], method="spearman")
# }
# simFuns <- c("spearman", "trunSpearman")
# simFunsDir <- c(F, F)
# tmp <- assessSig(simDat, simFuns, simFunsDir, nPerm=10, deg=49)

## test parallelization
# nCores <- parallel::detectCores()
# registerDoParallel(nCores)
# system.time({tmp <- foreach(iterators::icount(100), .combine='rbind') %dopar%
#   sapply(simFuns,
#          function(fun) {
#            do.call(fun, args=list(x=sample(dat[[1]][[1]]$x),
#                                   y=dat[[1]][[1]]$y,
#                                   deg=deg))
#          })})
# system.time({ tmp <- assessSig(dat[[1]][[1]], simFuns, simFunsDir, nPerm=10, deg=49) })
# stopImplicitCluster()
