# Some utility functions
nondiagpick <- function(x,i){
  # For spanning the \lambda - \lambda part of the Fisher Information Matrix
  n <- length(x)
  if (i == 1){
    return(c(x[1]/2,rep(0,n-i)))
  }
  return(c(rep(x[i],i-1),x[i]/2,rep(0,n-i)))
}

rcumsumr <- function(x){
  # for constructing \sum_{j >= i} something_j
  if(!is.vector(x)){
    x <- as.vector(x)
  }
  rev(cumsum(rev(x)))
}

lpNorm <- function(x,p){
  ss <- sum(x^p)
  (ss)^(1/p)
}