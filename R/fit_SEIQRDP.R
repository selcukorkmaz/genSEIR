fit_SEIQRDP <- function(Q, R, D, Npop, E0, I0, time, guess, ftol = 1e-6,
                        ptol = 1e-6, gtol = 1e-6, epsfcn = 0.001, factor = 100, maxfev = 1000,
                        maxiter = 100, nprint = 1, trace = TRUE, ...){

  #'Fit SEIQRDP function
  #'
  #'Fit SEIQRDP function parameters used in the SEIQRDP function, used to model
  #' the time-evolution of an epidemic outbreak.

  #' @param Q time histories of the active cases
  #' @param R time histories of the recovered cases
  #' @param D time histories of the deceased cases
  #' @param Npop total population of the country
  #' @param E0 initial number of exposed cases
  #' @param I0 initial number of predicted infectious cases
  #' @param time time vector
  #' @param guess initiail guess parameters
  #' @param ftol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
  #' @param ptol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
  #' @param gtol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
  #' @param epsfcn nls.lm.control object. Default is \code{0.001}
  #' @param factor nls.lm.control object. Default is \code{100}
  #' @param maxfev nls.lm.control object. Default is \code{1000}
  #' @param maxiter nls.lm.control object. Default is \code{100}
  #' @param nprint nls.lm.control object. Default is \code{1}
  #' @param trace set \code{TRUE} to trace iteration results
  #' @param ... further arguments
  #'
  #'@importFrom minpack.lm nls.lm nls.lm.control
  #'@export fit_SEIQRDP

  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}


  Q[Q<0] <- 0
  R[R<0] <- 0
  D[D<0] <- 0

  if (is.null(R)){
    warning(' No data available for "Recovered" ')
    input = cbind.data.frame(Q=as.numeric(Q), D=as.numeric(D)) #% In this aprticular case, Q is actually the number of active + recovered cases
  }else{
    input = cbind.data.frame(Q=as.numeric(Q),R=as.numeric(R),D=as.numeric(D))
  }

  if(is.null(time)){

    stop('Time should be a vector')
  }

  fs = 1/dt
  tTarget = as.numeric(round((time-time[1])*fs)/fs)
  t = seq(tTarget[1], tTarget[length(tTarget)], dt)

  if (!is.null(R)){
    getLambda = getLambdaFun(tTarget, Q, R, guess, ftol,
                             ptol, gtol, epsfcn, factor, maxfev,
                             maxiter, nprint, trace)
    guess = getLambda$guess
    lambdaFun = getLambda$lambdaFun
  }else{
    lambdaFun =  function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))}
  }

  getKappa = getKappaFun(tTarget, Q, D, guess, ftol,
                         ptol, gtol, epsfcn, factor, maxfev,
                         maxiter, nprint, trace)
  guess = getKappa$guess
  kappaFun = getKappa$kappaFun

  modelFun1 = SEIQRDP_for_fitting

  if (is.null(R)){
    kappaMax = guess[8:10]*1.05
    kappaMin = guess[8:10]*0.95
    lambdaMax = c(1, 1, 100)
    lambdaMin = c(0, 0, 0)

    if (kappaMax[3]<1e-1){
      kappaMax[3] = 100
      kappaMin[3] = 0
    }

  }else{
    kappaMax = guess[8:10]*2.0
    kappaMin = guess[8:10]/2.0
    lambdaMax = guess[5:7]*2.0
    lambdaMin = guess[5:7]/2.0

    if (kappaMax[3]<1e-1){
      kappaMax[3] = 20
      kappaMin[3] = 0
    }

    if (lambdaMax[3]<1e-1){
      lambdaMax[3] = 20
      lambdaMin[3] = 0
    }
  }

  ub = c(1, 5, 1, 1, lambdaMax, kappaMax)
  lb = c(0, 0, 0, 0, lambdaMin, kappaMin)

  df = cbind.data.frame(input)
  colnames(df) = c("Q1","R1","D1")

  ctrl = nls.lm.control(ftol = ftol,
                        ptol = ptol, gtol = gtol, diag = list(), epsfcn = epsfcn,
                        factor = factor, maxfev = maxfev, maxiter = maxiter, nprint = nprint)

  func <- function(p, obs, t11,  t00) {obs - SEIQRDP_for_fitting(p, t11,  t00)}

  nls.out <- nls.lm(par=guess, fn = func,  t11 = t,
                    t00=tTarget, obs =data.matrix(df),
                    lower = lb, upper = ub, control = ctrl)

  Coeff = nls.out$par

  alpha1 = abs(Coeff[1])
  beta1 = abs(Coeff[2])
  gamma1 = abs(Coeff[3])
  delta1 = abs(Coeff[4])
  Lambda1 = abs(Coeff[5:7])
  Kappa1 = abs(Coeff[8:10])

  return(list(alpha1=alpha1,beta1=beta1,gamma1=gamma1,delta1=delta1,
              Lambda1=Lambda1,Kappa1=Kappa1,lambdaFun=lambdaFun,kappaFun=kappaFun))

}
