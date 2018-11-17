
assessSig <- function(expDat, simFuns, nPerm, df.x) {
  ###  computes similarity scores for nPerm permutations
  # expDat <- expression data: data.frame with x = signature1 (disease) and y=signature2 (drug)
  # simFuns: vector of similarity functions
  # nPerm: number of permutations
  permResults <- sapply(1:nPerm,
                        function(x) {
                          sapply(simFuns,
                                 function(fun) {
                                   do.call(fun, args=list(x=sample(expDat$x),
                                                          y=expDat$y,
                                                          df.x=df.x))
                                 })
                        })
  names(permResults) <- simFuns
  return(t(permResults))
}

## example
simDat <- data.frame(x=rt(20000, df=49)+log(5), y=rt(20000, df=49))
spearman <- function(x, y, df.x) {
  cor(x, y, method="spearman")
}
trunSpearman <- function(x, y, df.x) {
  x.pvalues <- pt(x, df=df.x)
  x.pvalues[x.pvalues>0.5] <- 1-x.pvalues[x.pvalues>0.5]
  x.pvalues <- 2*x.pvalues
  use <- p.adjust(x.pvalues, method="fdr")<0.05
  cor(x[use], y[use], method="spearman")
}
simFuns <- c("spearman", "trunSpearman")
tmp <- assessSig(simDat, simFuns, 10, df.x=49)
