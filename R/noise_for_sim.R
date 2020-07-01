noise <- function(X,s){
  out <- rnorm(n = length(X),mean = X, sd = s)
  return(out)
}
