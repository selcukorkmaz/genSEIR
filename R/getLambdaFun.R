#' Estimate Death Rate
#'
#' This function provides a first estimate of the  recovery rate, to faciliate
#' convergence of the main algorithm.

#' @param tTarget target time vector
#' @param Q target time-histories of the quarantined cases
#' @param R target time-histories of the recovered cases
#' @param guess Initiail guess parameters for kappa
#' @param ftol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param ptol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param gtol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param epsfcn nls.lm.control object. Default is \code{0.001}
#' @param factor nls.lm.control object. Default is \code{100}
#' @param maxfev nls.lm.control object. Default is \code{1000}
#' @param maxiter nls.lm.control object. Default is \code{100}
#' @param nprint nls.lm.control object. Default is \code{1}
#' @param trace set \code{TRUE} to trace iteration results
#'
#' @importFrom nlsr nlxb
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return vector of estimation and optimization function for the recovery rate
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

getLambdaFun <- function (tTarget, Q, R, guess, ftol,
                          ptol, gtol, epsfcn, factor, maxfev,
                          maxiter, nprint, trace){

  if (max(R)<20){
    lambdaFun =  function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))}
  }else{

  myFun1 = function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))};
  myFun2 = function(a,t) {a[1] + exp(-a[2]*(t+a[3]))};

  rate = diff(as.numeric(R))/median(diff(tTarget))/as.numeric(Q)[2:length(as.numeric(Q))]
  x = tTarget[2:length(tTarget)]

  rate[abs(rate)>1 | abs(rate)==0]=NA

  df = cbind.data.frame(tk = x[!is.na(rate)], z = rate[!is.na(rate)])


  ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)

  model1 <- nlxb(z ~ a1 / (1+exp(-a2*(tk-a3))),
                  start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
                 control = ctrl, lower = c(0,0,0), upper = c(1,1,100))

  coeff1 = model1$coefficients
  r1 = model1$ssquares

  ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)


  model2 <- nlxb(z ~ a1 + exp(-a2*(tk+a3)),
                  start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
                 control = ctrl, lower = c(0,0,0), upper = c(1,1,100))
  coeff2 = model2$coefficients
  r2 = model2$ssquares

  if (r1<r2 || coeff2[1]>0.99 || coeff2[2]>4.9){
    lambdaGuess = coeff1
    lambdaFun = myFun1
  }else{
    lambdaGuess = coeff2
    lambdaFun = myFun2
  }

  guess[5:7] = lambdaGuess

  }
  return(list(guess=guess, lambdaFun = lambdaFun))

}
