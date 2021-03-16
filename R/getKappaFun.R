kappaFun = function(a,t) {a[1] / (exp(a[2]*(t-a[3])))}

getKappaFun <- function(tTarget, Q, D, guess, ftol,
                        ptol, gtol, epsfcn, factor, maxfev,
                        maxiter, nprint, trace){

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

  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

    if (max(D)<10){
    kappaFun = function(a,t) {a[1] / (exp(a[2]*(t-a[3])))}

    }else{

    #   try
    # opt=optimset('TolX',1e-6,'TolFun',1e-6,'Display','off');

    # %                myFun1 = @(a,t) a(1).*exp(-a(2)*(t+(a(3))));

    myFun1 = function(a,t) {a[1] / (exp(a[2]*(t-a[3])) + exp(-a[2]*(t-a[3])))}
    myFun2 = function(a,t) {a[1]*exp(-(a[2]*(t-a[3]))^2)}
    myFun3 = function(a,t) {a[1] + exp(-a[2]*(t+a[3]))}

    rate = as.numeric((diff(as.numeric(D))/median(diff(tTarget)))/Q[2:length(Q)])
    x = tTarget[2:length(tTarget)]

    # % A death rate larger than 3 is abnormally high. It is not
    # % used for the fitting.
    rate[abs(rate)>3]=NA

    # % Remove death rate = 0 if the majority number of death is not
    # % zero
    if (length(rate[rate==0])/length(rate) < 0.5){
      rate[abs(rate)==0] = NA
    }


    df = cbind.data.frame(tk = x[!is.na(rate)], z = rate[!is.na(rate)])

    # lsqCurveFit1 = pracma::lsqcurvefit(myFun1, guess[8:10], x[!is.na(rate)], rate[!is.na(rate)])
    # coeff1 = lsqCurveFit1$x
    # r1 = lsqCurveFit1$ssq
    #
    ctrl = nls.lm.control(ftol = ftol,
                          ptol = ptol, gtol = gtol, diag = list(), epsfcn = epsfcn,
                          factor = factor, maxfev = maxfev, maxiter = maxiter, nprint = nprint)

    # ctrlNLS <- nls.control(maxiter = 2500, tol = 1e-05, minFactor = 1/1024,
    #             printEval = T, warnOnly = T)
    #
    #
    # lsqCurveFit1 <- nlsLM(z ~ (a1 / (exp(a2*(tk-a3)) + exp(-a2*(tk-a3)))),
    #                 start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace= trace,
    #                 control=ctrl, lower = c(0,0,0), upper = c(1,1,100))
    #
    # s = summary(lsqCurveFit1)
    # coeff1 = round(s$parameters[,1],4)
    # r1 = sum((s$residuals)^2)
    # coeff1
    # r1

    #
    ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)
    lsqCurveFit1 <- nlxb(z ~ (a1 / (exp(a2*(tk-a3)) + exp(-a2*(tk-a3)))),
                          start=c(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace= trace,
                          control = ctrl, lower = c(0,0,0), upper = c(1,1,100))


    coeff1 = lsqCurveFit1$coefficients
    r1 = lsqCurveFit1$ssquares


#
#     lsqCurveFit2 <- nlsLM(z ~ a1*exp(-(a2*(tk-a3))^2),
#                           start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
#                           control=ctrl, lower = c(0,0,0), upper = c(1,1,100))
#     s2 = summary(lsqCurveFit2)
#     coeff2 = s2$parameters[,1]
#     r2 = sum((s2$residuals)^2)

    ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)

    lsqCurveFit2 <- nlxb(z ~ a1*exp(-(a2*(tk-a3))^2),
                          start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
                          control = ctrl, lower = c(0,0,0), upper = c(1,1,100))
    coeff2 = lsqCurveFit2$coefficients
    r2 = lsqCurveFit2$ssquares


    # lsqCurveFit3 <- nlsLM(z ~ a1 + exp(-a2*(tk+a3)),
    #                       start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
    #                       control=ctrl, lower = c(0,0,0), upper = c(1,1,100))
    # s3 = summary(lsqCurveFit3)
    # coeff3 = s3$parameters[,1]
    # r3 = sum((s3$residuals)^2)

    lsqCurveFit3 <- nlxb(z ~ a1 + exp(-a2*(tk+a3)),
                          start=list(a1=guess[8], a2=guess[9], a3=guess[10]), data=df, trace=trace,
                         control = ctrl, lower = c(0,0,0), upper = c(1,1,100))

    coeff3 = lsqCurveFit3$coefficients
    r3 = lsqCurveFit3$ssquares

    # %
    # %                     figure;plot(x,rate,x,myFun1(coeff1,x),'r',x,myFun2(coeff2,x),'g',x,myFun3(coeff3,x),'b')

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

    guess[8:10] = kappaGuess # update guess

    }

    return(list(guess=guess,kappaFun=kappaFun))

}
