#' Compute the matrix A
#'
#' This function computes the matrix A that is found in: dY/dt = A*Y + F

#' @param alpha protection rate
#' @param gamma inverse of the average latent time
#' @param delta rate of people entering in quarantine
#' @param lambda cure rate
#' @param kappa mortality rate
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#'@return The matrix A that is found in: dY/dt = A*Y + F
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

getA <- function(alpha, gamma, delta, lambda, kappa){

  A = matrix(0,7,7)
  # S
  A[1,1] = -alpha
  # E
  A[2,2] = -gamma
  # I
  A[3,2:3] = c(gamma,-delta)
  # Q
  A[4,3:4] = c(delta,-kappa-lambda)
  # R
  A[5,4] = lambda
  # D
  A[6,4] = kappa
  # P
  A[7,1] = alpha

  return(A)
}
