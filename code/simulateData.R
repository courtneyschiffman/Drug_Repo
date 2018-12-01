
simulateData <- function(p1, p2, p3, p4, nGenes=20000, mu1=log(8), mu2=log(8), df1=49, df2=49,
                         frac1=0.5, frac2=0.5) {
  ### simulates pair of vectors of t-statistics representing 2 gene expression signatures
  ### returns: data.frame with x = signature1 (disease) and y=signature2 (drug)
  # p1: fraction of signature1 that is differentially expressed
  # p2: fraction of signature2 that is differentially expressed
  # p3: fraction of DE in sig1 that is also DE in same direction in sig2
  # p4: fraction of DE in sig1 that is also DE but in different direction in sig2
  # n: number of genes
  # mu1: mean of alternative hypothesis for signature1
  # mu2: mean of alternative hypothesis for signature2
  # df1: degrees of freedom for signature1
  # df2: digrees of freedom for signature2
  # frac1: fraction of upregulated DE genes in signature 1
  # frac2: fraction of upregulated DE genes in signature 2

  ## check p3+p4≤1
  if(p3+p4>1) {
    print("invalid p3 and p4 values")
    return(NULL)
  }

  ## signature 1
  nDE.1 <- round(nGenes*p1)
  nDE.1.up <- round(nDE.1*frac1)

  ## signature 2
  nDE.2.up <- round(nGenes*p2*frac2)
  nDE.2.up.same <- round(nDE.2.up*p3)
  nDE.2.up.diff <- round(nDE.2.up*p4)
  nDE.2.down <- round(nGenes*p2*(1-frac2))
  nDE.2.down.same <- round(nDE.2.down*p3)
  nDE.2.down.diff <- round(nDE.2.down*p4)

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

  sig2 <- c(rt(nDE.2.up.same, df2)+mu2, # upregulated in sig1 and sig2
            rt(nDE.2.down.diff, df2)-mu2, # upregulated in sig1, downregulated in sig2
            rt(nDE.1.up-nDE.2.up.same-nDE.2.up.diff, df2), # upregulated in only sig1
            rt(nDE.2.up.diff, df2)+mu2, # downregulated in sig1, upregulated in sig2
            rt(nDE.2.down.same, df2)-mu2, # downregulated in sig1 and sig2
            rt((nDE.1-nDE.1.up)-nDE.2.up.diff-nDE.2.down.same, df2), # downregulated in only sig1
            rt(nDE.2.up-nDE.2.up.same-nDE.2.up.diff, df2)+mu2, # upregulated in only sig2
            rt(nDE.2.down-nDE.2.down.same-nDE.2.down.diff, df2)-mu2, # downregulated in only sig2
            rt((nGenes-nDE.1)-(nDE.2.up-nDE.2.up.same-nDE.2.up.diff)-(nDE.2.down-nDE.2.down.same-nDE.2.down.diff), df2))
  return(data.frame(x=sig1, y=sig2))
}

# # example 1
# tmp <- simulateData(p1=0.05, p2=0.05, p3=0.1, p4=0.8)
# indices <- as.list(rep(NA, 6))
# names(indices) <- apply(expand.grid(c("sig1", "sig2"), c("up", "down", "null")),1,paste0, collapse=" ")
# indices[[1]] <- c(1:500)
# indices[[2]] <- c(1:50, 501:900, 1001:1050)
# indices[[3]] <- c(501:1000)
# indices[[4]] <- c(51:450, 901:950, 1051:1100)
# indices[[5]] <- c(1001:20000)
# indices[[6]] <- c(451:500, 951:1000, 1101:20000)
# tmpSpearman <- sapply(indices, function(i) { cor(x=tmp$x[i], y=tmp$y[i], method="spearman") })
# par(mfrow=c(3,2))
# for(i in 1:length(indices)) {
#   plot(tmp$x[indices[[i]]], tmp$y[indices[[i]]], pch=19, xlab="", ylab="",
#        main=paste(names(indices)[i], "; spearman =", round(tmpSpearman[i],2)))
#   abline(a=0, b=1, col="red")
#   abline(h=0, col="blue")
#   abline(v=0, col="blue")
# }
#
# # example 2
# tmp <- simulateData(p1=0.05, p2=0.05, p3=0.2, p4=0.4)
# indices <- as.list(rep(NA, 6))
# names(indices) <- apply(expand.grid(c("sig1", "sig2"), c("up", "down", "null")),1,paste0, collapse=" ")
# indices[[1]] <- c(1:500)
# indices[[2]] <- c(1:50, 501:900, 1001:1050)
# indices[[3]] <- c(501:1000)
# indices[[4]] <- c(51:450, 901:950, 1051:1100)
# indices[[5]] <- c(1001:20000)
# indices[[6]] <- c(451:500, 951:1000, 1101:20000)
# tmpSpearman <- sapply(indices, function(i) { cor(x=tmp$x[i], y=tmp$y[i], method="spearman") })
# par(mfrow=c(3,2))
# for(i in 1:length(indices)) {
#   plot(tmp$x[indices[[i]]], tmp$y[indices[[i]]], pch=19, xlab="", ylab="",
#        main=paste(names(indices)[i], "; spearman =", round(tmpSpearman[i],2)))
#   abline(a=0, b=1, col="red")
#   abline(h=0, col="blue")
#   abline(v=0, col="blue")
# }
