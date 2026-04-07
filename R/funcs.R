inv_logistic_transf <- function(x) {
  return(2 * (exp(x)/(1+exp(x))) - 1)
}

logistic_transf<- function(x) {
  log((1+x)/(1-x))
}

rollmean_4yr <- function(x) {
  stats::filter(x, rep(1/4, 4), sides = 1)
}
