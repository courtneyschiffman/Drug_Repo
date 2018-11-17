
simulateData <- function(p1, p2, p3, nGenes=20000, mu1=log(5), mu2=log(5), df1=49, df2=49) {
  ### simulates pair of vectors of t-statistics representing 2 gene expression signatures
  ### returns: data.frame with x = signature1 (disease) and y=signature2 (drug)
  # p1: fraction of signature1 that is differentially expressed
  # p2: fraction of signature2 that is differentially expressed
  # p3: fraction of DE signature2 that is also DE in signature1
  # n: number of genes
  # mu1: mean of alternative hypothesis for signature1
  # mu2: mean of alternative hypothesis for signature2
  # df1: degrees of freedom for signature1
  # df2: digrees of freedom for signature2
  nDE.1 <- round(nGenes*p1)
  nDE.2 <- round(nGenes*p2)
  nDE.1and2 <- round(nDE.2*p3)
  sig1 <- c(rt(nDE.1, df1)+mu1, # DE genes
            rt(nGenes-nDE.1, df1)) # non-DE genes
  sig2 <- c(rt(nDE.1and2, df2)+mu2, # genes that are DE in sig1 and sig2
            rt(nDE.1-nDE.1and2, df2), # genes that are only DE in sig1
            rt(nDE.2-nDE.1and2, df2) + mu2, # genes that are only DE in sig2
            rt(nGenes-nDE.1-nDE.2+nDE.1and2, df2)) # non-DE genes
  return(data.frame(x=sig1, y=sig2))
}

## example
# tmp <- simulateData(p1=0.05, p2=0.05, p3=0)
# par(mfrow=c(2,1))
# plot(density(tmp$sig1), xlim=c(min(tmp), max(tmp)), main="signature1")
# abline(v=mean(tmp$sig1), col="red")
# plot(density(tmp$sig2), xlim=c(min(tmp), max(tmp)), main="signature2")
# abline(v=mean(tmp$sig2), col="red")
#
