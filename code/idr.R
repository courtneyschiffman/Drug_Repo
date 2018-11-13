
### SOME EXAMPLES

x <- c(rt(4000,df=49)+4,rt(4000,df=49)-4,rt(8000,df=49))
y <- c(-x[1:8000],rt(8000,df=49))

x <- c(rt(16000,df=49))
y <- c(rt(16000,df=49))


## x and y are your vectors of t-statistics

## Just made up a degrees of freedom here for the t-distributions to get the p-values
## I'm using correlation among the reproducible p-values as the statistic
## in the permutation tests look for any pairs that have a correlation higher than the obersved correlation
## correlation calculation is really robust to choices of mu and signma


IDR.func <- function(x,y,deg=49,mu=6,sigma=1.5,rho=.6,p=.2){
  which.tail <- c(TRUE,FALSE)
  pval.func <- function(i){
    m <- which.min(c(pt(x[i],df=deg,lower.tail = TRUE),pt(x[i],df=deg,lower.tail = FALSE)))
    x.pval <- pt(x[i],df=deg,lower.tail = which.tail[m])
    y.pval <- pt(y[i],df=deg,lower.tail = which.tail[-m])
    return(c(x.pval,y.pval))
  }
  
  new.pvals <- sapply(1:length(x),pval.func)
  new.pvals <- t(-log(new.pvals))
  my.idr <- est.IDR(new.pvals,mu=mu,sigma=sigma,rho=rho,p=p)$para$rho
  return(my.idr)
}

IDR.func(x,y)

IDR.func(sample(x),y)

plot(rank(new.pvals[,1]),rank(new.pvals[,2]))
