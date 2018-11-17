

spearman <- function(x, y, df.x) {
  ### returns spearman rank correlation for 2 vectors of t-statistics
  # x: vector of t-statistics (signature1)
  # y: vector of t-statistics (signature2)
  # df.x: degrees of freedom for signature1
  cor(x, y, method="spearman")
}

trunSpearman <- function(x, y, df.x) {
  ### returns spearman rank correlation only on DE genes
  # x: vector of t-statistics (signature 1)
  # y: vector of t-statistics (signature 2)
  # df.x: degrees of freedom for signature1

  ## compute pvalues for 2-sided tests
  x.pvalues <- pt(x, df=df.x)
  x.pvalues[x.pvalues>0.5] <- 1-x.pvalues[x.pvalues>0.5]
  x.pvalues <- 2*x.pvalues
  ## compute spearman rank correlation on only DE genes
  ## use FDR for multiple hypothesis correction
  use <- p.adjust(x.pvalues, method="fdr")<0.05
  cor(x[use], y[use], method="spearman")
}
