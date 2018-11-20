## x and y are vectors of t statistics.

## x is a vector of t statistics for the disease and likewise y is for the drug

x <- c(rt(4000,df=49)+4,rt(4000,df=49)-4,rt(8000,df=49))
y <- c(x[1:8000],rt(8000,df=49))


gsea <- function(x,y,deg){
  names(y) <- as.character(1:length(y))
  names(x) <- as.character(1:length(x))
  y <- rank(-y)
  xpvals <- 2*(1-pt(abs(x),df=deg,lower.tail = T))
  xqvals <- qvalue(xpvals)$qvalues
  xup <- sort(x[x>=0 & xqvals <= .05],decreasing = T)
  xdown <- sort(x[x<0 & xqvals <= .05],decreasing = T)
  vup <- y[names(y)%in%names(xup)]
  vdown <- y[names(y)%in%names(xdown)]
  vup <- vup[order(match(names(vup),names(xup)))]
  vdown <- vdown[order(match(names(vdown),names(xdown)))]
  
  aup <- max(sapply(1:length(xup),function(j) (j/length(xup)-vup[j]/length(y))))
  bup <- max(sapply(1:length(xup),function(j) (vup[j]/length(y)-(j-1)/length(xup))))
  adown <- max(sapply(1:length(xdown),function(j) (j/length(xdown)-vdown[j]/length(y))))
  bdown <- max(sapply(1:length(xdown),function(j) (vdown[j]/length(y)-(j-1)/length(xdown))))
  if(aup>bup){
    esup <- aup
  }else{
    esup <- -bup
  }
  if(adown>bdown){
    esdown <- adown
  }else{
    esdown <- -bdown
  }
  if(sign(esup)==sign(esdown)){
    dds <- 0
  }else{
    dds <- esup-esdown
  }
  return(dds)
}

x <- c(rt(4000,df=49)+4,rt(4000,df=49)-4,rt(8000,df=49))
y <- c(x[1:8000],rt(8000,df=49))

gsea(x,sample(y),49)

x <- c(rt(4000,df=49)+4,rt(4000,df=49)-4,rt(8000,df=49))
y <- c(-x[1:8000],rt(8000,df=49))
gsea(x,y,49)

y <- rt(20000,df=49) + 4
x <- rt(20000,df=49) + 4

gsea(x,y,49)

y <- rt(20000,df=49) + 4
x <- -y

gsea(x,y,49)


