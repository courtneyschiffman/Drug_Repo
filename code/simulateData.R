
simulateData <- function(p1, p2, p3, nGenes=20000, mu1=log(8), mu2=log(8), df1=49, df2=49,
                         frac1=0.5, frac2=0.5, frac3=0.5) {
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
  # frac1: fraction of upregulated DE genes in signature 1
  # frac2: fraction of upregulated DE genes in signature 2
  # frac3: fraction of upregulated DE genes in signatures 1 and 2
  nDE.1 <- round(nGenes*p1)
  nDE.1.up <- round(nDE.1*frac1)
  nDE.2 <- round(nGenes*p2)
  nDE.2.up <- round(nDE.2*frac2)
  nDE.1and2 <- round(nDE.2*p3)
  nDE.1and2.up <- round(nDE.1and2*frac3)
  checkDE <- T
  while(checkDE) {
    sig1 <- c(rt(nDE.1.up, df1)+mu1, # upregulated in sig1
              rt(nDE.1-nDE.1.up, df1)-mu2, # downregualted in sig1
              rt(nGenes-nDE.1, df1)) # non-DE genes
    qvals <- p.adjust(2*pt(abs(sig1), df=df1, lower.tail=F), method="fdr")
    sig1.nUp <- sum(sig1>0 & qvals<=0.05)
    sig1.nDown <- sum(sig1<0 & qvals<=0.05)
    checkDE <- !(sig1.nUp>0 & sig1.nDown>0)
  }
  sig2 <- c(rt(nDE.1and2.up, df2)+mu2, # upregulated in sig1 and sig2
            rt(nDE.1.up-nDE.1and2.up, df2), # upregulated in only sig1
            rt(nDE.1and2-nDE.1and2.up, df2)-mu2, # downregulated in sig1 and sig2
            rt(nDE.1-nDE.1and2-(nDE.1.up-nDE.1and2.up), df2), # downregulated DE in only sig2
            rt(nDE.2.up-nDE.1and2.up, df2)+mu2, # upregulated in only sig2
            rt((nDE.2-nDE.2.up)-(nDE.1and2-nDE.1and2.up), df2)-mu2, # downregulated in only sig2
            rt((nGenes-nDE.1)-(nDE.2-nDE.1and2), df2)) # non-DE genes
  return(data.frame(x=sig1, y=sig2))
}

## example
tmp <- simulateData(p1=0.2, p2=0.05, p3=0)
pvalues <- p.adjust(2*pt(abs(tmp$x), df=49, lower.tail=F), method="fdr")
sum(tmp$x<0&pvalues<=0.05)
sum(tmp$x>0&pvalues<=0.05)
# par(mfrow=c(2,1))
# plot(density(tmp$x), xlim=c(min(tmp), max(tmp)), main="signature1")
# abline(v=mean(tmp$x), col="red")
# plot(density(tmp$y), xlim=c(min(tmp), max(tmp)), main="signature2")
# abline(v=mean(tmp$y), col="red")

