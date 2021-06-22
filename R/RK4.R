#' Runge-Kutta 4th Order Method to Solve Differential Equation
#'
#' @param Y time vector
#' @param A the matrix A that is found in: dY/dt = A*Y + F
#' @param K the zero matrix for the seven states
#' @param dt the time step. This oversamples time to ensure that the algorithm converges
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return ordinary differential equation result as a vector
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

RK4 <- function(Y, A, K, dt){

  k_1 = modelFun(Y, A, K)
  k_2 = modelFun(Y+0.5*dt*k_1, A, K)
  k_3 = modelFun(Y+0.5*dt*k_2, A, K)
  k_4 = modelFun(Y+k_3*dt, A, K)
  Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt
  return(Y)
}
