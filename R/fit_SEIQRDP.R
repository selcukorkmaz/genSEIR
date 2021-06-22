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
#' @param time a time vector
#' @param dt the time step. This oversamples time to ensure that the algorithm converges
#' @param guess initial guess parameters
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
#'
#'@export fit_SEIQRDP
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return a list of optimized parameters
#'
#' @examples
#'\donttest{
#'start = "01/01/21"
#'finish = "04/01/21"
#'country = "Italy"
#'dt = 1
#'f=30
#'
#'covidData = getDataCOVID(start = start, finish = finish, country = country)
#'Recovered = covidData$tableRecovered
#'Deaths = covidData$tableDeaths
#'Confirmed = covidData$tableConfirmed
#'
#'if(nrow(Recovered) == 1){
#'   name = Recovered$CountryRegion
#'}else{
#'    name = paste0(Recovered$ProvinceState, " (",Recovered$CountryRegion,")")
#'}
#'
#'   recovered = Recovered[ ,5:ncol(covidData$tableRecovered)]
#'   deaths = Deaths[ ,5:ncol(covidData$tableDeaths)]
#'   confirmed = Confirmed[ ,5:ncol(covidData$tableConfirmed)]
#'
#'   Npop = 60000000
#'
#'   alpha_guess = 0.05
#'   beta_guess = 0.8
#'   LT_guess = 7
#'   Q_guess = 0.8
#'   lambda_guess = c(0.01,0.001,10)
#'   kappa_guess = c(0.001,0.001,10)
#'
#'   guess = c(alpha_guess,
#'             beta_guess,
#'             1/LT_guess,
#'             Q_guess,
#'             lambda_guess,
#'             kappa_guess)
#'
#'  Q0 = confirmed[1]-recovered[1]-deaths[1]
#'  I0 = 0.3*Q0
#'  E0 = 0.3*Q0
#'  R0 = recovered[1]
#'  D0 = deaths[1]
#'
#'  Active = confirmed-recovered-deaths
#'  Active[Active<0] <- 0
#'
#'  Q=Active
#'  R=recovered
#'  D = deaths
#'
#'  time = seq(as.Date(start, format = "%m/%d/%y"), as.Date(finish, format = "%m/%d/%y"), by = "1 day")
#'
#'  params = fit_SEIQRDP(Q = Active, R = recovered, D = deaths, Npop = Npop, E0 = E0, I0 = I0,
#'                         time = time, dt = dt, guess = guess, ftol = 1e-6, ptol = 1e-6, gtol = 1e-6,
#'                         epsfcn = 0.001, factor = 100, maxfev = 1000,maxiter = 100, nprint = 1,
#'                         trace = TRUE)
#'}
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

fit_SEIQRDP <- function(Q, R, D, Npop, E0, I0, time, dt = 1/24, guess, ftol = 1e-6,
                        ptol = 1e-6, gtol = 1e-6, epsfcn = 0.001, factor = 100, maxfev = 1000,
                        maxiter = 100, nprint = 1, trace = TRUE, ...){

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
