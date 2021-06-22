#'Model function
#'
#' @param Y time vector
#' @param A the matrix A that is found in: dY/dt = A*Y + F
#' @param K the zero matrix for the seven states
#'
#' @return No return value, called for side effects

modelFun = function(Y, A, K) {A %*% Y + K}
