#' Anonymous function approximating the death rate
#'
#' @param a parameter vector
#' @param t time vector
#'
#'@return No return value, called for side effects

kappaFun = function(a, t) {a[1] / (exp(a[2]*(t-a[3])))}
