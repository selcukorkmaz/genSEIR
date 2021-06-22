#' Estimate Death Rate
#'
#' This function provides a first estimate of the  death rate, to faciliate
#' convergence of the main algorithm.

#' @param tTarget time vector
#' @param Q target time-histories of the quarantined cases
#' @param D target time-histories of the dead cases
#' @param guess Initiail guess parameters for kappa
#' @param ftol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param ptol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param gtol nls.lm.control object. non-negative numeric. Default is \code{1e-6}
#' @param epsfcn nls.lm.control object. Default is \code{0.001}
#' @param factor nls.lm.control object. Default is \code{100}
#' @param maxfev nls.lm.control object. Default is \code{1000}
#' @param maxiter nls.lm.control object. Default is \code{100}
#' @param nprint nls.lm.control object. Default is \code{1}
#' @param trace Set \code{TRUE} to trace iteration results
#'
#' @importFrom nlsr nlxb
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return vector of estimation and optimization function for the death rate
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}


getKappaFun <- function(tTarget, Q, D, guess, ftol,
                        ptol, gtol, epsfcn, factor, maxfev,
                        maxiter, nprint, trace){

    if (max(D)<10){
    kappaFun = function(a,t) {a[1] / (exp(a[2]*(t-a[3])))}

    }else{

    myFun1 = function(a,t) {a[1] / (exp(a[2]*(t-a[3])) + exp(-a[2]*(t-a[3])))}
    myFun2 = function(a,t) {a[1]*exp(-(a[2]*(t-a[3]))^2)}
    myFun3 = function(a,t) {a[1] + exp(-a[2]*(t+a[3]))}

    rate = as.numeric((diff(as.numeric(D))/median(diff(tTarget)))/Q[2:length(Q)])
    x = tTarget[2:length(tTarget)]

    rate[abs(rate)>3]=NA

    if (length(rate[rate==0])/length(rate) < 0.5){
      rate[abs(rate)==0] = NA
    }

    df = cbind.data.frame(tk = x[!is.na(rate)], z = rate[!is.na(rate)])

    ctrl = nls.lm.control(ftol = ftol,
                          ptol = ptol, gtol = gtol, diag = list(), epsfcn = epsfcn,
                          factor = factor, maxfev = maxfev, maxiter = maxiter, nprint = nprint)


    ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)
    lsqCurveFit1 <- nlxb(z ~ (a1 / (exp(a2*(tk-a3)) + exp(-a2*(tk-a3)))),
                          start=c(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace= trace,
                          control = ctrl, lower = c(0,0,0), upper = c(1,1,100))


    coeff1 = lsqCurveFit1$coefficients
    r1 = lsqCurveFit1$ssquares


    ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)

    lsqCurveFit2 <- nlxb(z ~ a1*exp(-(a2*(tk-a3))^2),
                          start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
                          control = ctrl, lower = c(0,0,0), upper = c(1,1,100))
    coeff2 = lsqCurveFit2$coefficients
    r2 = lsqCurveFit2$ssquares


    lsqCurveFit3 <- nlxb(z ~ a1 + exp(-a2*(tk+a3)),
                          start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
                         control = ctrl, lower = c(0,0,0), upper = c(1,1,100))

    coeff3 = lsqCurveFit3$coefficients
    r3 = lsqCurveFit3$ssquares

    minR = min(r1,r2,r3)
    if (r1==minR){
      kappaGuess = coeff1
      kappaFun = myFun1
    }else if (r2==minR){
      kappaGuess = coeff2
      kappaFun = myFun2
    }else if (r3==minR){
      kappaFun = myFun3
      kappaGuess = coeff3
    }

    guess[8:10] = kappaGuess

    }

    return(list(guess=guess, kappaFun=kappaFun))

}
