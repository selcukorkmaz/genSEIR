#' Anonymous function approximating the recovery rate
#'
#' @param a parameter vector
#' @param t time vector
#'
#' @return No return value, called for side effects

lambdaFun =  function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))}
