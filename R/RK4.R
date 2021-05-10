RK4 <- function(Y, A, K, dt){

  #' Runge-Kutta of order 4
  #'

  #' @param Y time vector
  #' @param A target time-histories of the quarantined cases
  #' @param K target time-histories of the recovered cases
  #' @param dt Initiail guess parameters for kappa


  #' @author Selcuk KORKMAZ, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

  k_1 = modelFun(Y, A, K)
  k_2 = modelFun(Y+0.5*dt*k_1, A, K)
  k_3 = modelFun(Y+0.5*dt*k_2, A, K)
  k_4 = modelFun(Y+k_3*dt, A, K)
  Y = Y + (1/6)*(k_1+2*k_2+2*k_3+k_4)*dt
  return(Y)
}
