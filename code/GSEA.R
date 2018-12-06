## x and y are vectors of t statistics.

## x is a vector of t statistics for the disease and likewise y is for the drug

gsea <- function(x,y,deg){
  names(y) <- as.character(1:length(y))
  names(x) <- as.character(1:length(x))
  y <- rank(-y)
  xpvals <- 2*pt(abs(x), df=deg, lower.tail=F)
  xqvals <- p.adjust(xpvals, method="fdr")
  xup <- sort(x[x>=0 & xqvals <= .05],decreasing = T)
  xdown <- sort(x[x<0 & xqvals <= .05],decreasing = T)
  vup <- y[names(y)%in%names(xup)]
  vdown <- y[names(y)%in%names(xdown)]
  vup <- sort(vup)
  vdown <- sort(vdown)
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


