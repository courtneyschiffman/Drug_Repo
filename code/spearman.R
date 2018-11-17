

spearman <- function(x, y, deg) {
  ### returns spearman rank correlation for 2 vectors of t-statistics
  # x: vector of t-statistics (signature1)
  # y: vector of t-statistics (signature2)
  # df.x: degrees of freedom for signature1
  cor(x, y, method="spearman")
}

trunSpearman <- function(x, y, deg) {
  ### returns spearman rank correlation only on DE genes
  # x: vector of t-statistics (signature 1)
  # y: vector of t-statistics (signature 2)
  # df.x: degrees of freedom for signature1
  x.pvalues <- 2*pt(abs(x), df=deg, lower.tail=F)
  use <- p.adjust(x.pvalues, method="fdr")<0.05
  cor(x[use], y[use], method="spearman")
}
