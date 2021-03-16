lambdaFun =  function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))}

getLambdaFun <- function (tTarget, Q, R, guess, ftol,
                          ptol, gtol, epsfcn, factor, maxfev,
                          maxiter, nprint, trace){


  #' Estimate Death Rate
  #'
  #' This function provides a first estimate of the  death rate, to faciliate
  #' convergence of the main algorithm.

  #' @param tTarget time vector
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
  #' @param trace Set \code{TRUE} to trace iteration results
  #'
  #' @importFrom nlsr nlxb

  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

  if (max(R)<20){
    lambdaFun =  function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))}
  }else{

    # try

  # opt=optimset('TolX',1e-6,'TolFun',1e-6,'Display','off');

  # % Two empirical functions are evaluated
  myFun1 = function(a,t) {a[1] / (1+exp(-a[2]*(t-a[3])))};
  myFun2 = function(a,t) {a[1] + exp(-a[2]*(t+a[3]))};

  # % Compute the recovery rate from the data (noisy data)
  rate = diff(as.numeric(R))/median(diff(tTarget))/as.numeric(Q)[2:length(as.numeric(Q))]
  x = tTarget[2:length(tTarget)]

  # % A daily rate larger than one is abnormally high. It is not
  # % used for the fitting. A daily recovered rate of zero is
  # % either abnormally low or reflects an insufficient number
  # % of recovered cases. It is not used either for the fitting
  rate[abs(rate)>1 | abs(rate)==0]=NA

  df = cbind.data.frame(tk = x[!is.na(rate)], z = rate[!is.na(rate)])

  # coeff1 = pracma::lsqcurvefit(myFun1, guess[5:7], x[!is.na(rate)], rate[!is.na(rate)])$x
  #r1 = pracma::lsqcurvefit(myFun1, guess[5:7], x[!is.na(rate)], rate[!is.na(rate)])$ssq


  # ctrl = nls.lm.control(ftol = ftol,
  #                       ptol = ptol, gtol = gtol, diag = list(), epsfcn = epsfcn,
  #                       factor = factor, maxfev = maxfev, maxiter = maxiter, nprint = nprint)

  # model1 <- nlsLM(z ~ a1 / (1+exp(-a2*(tk-a3))),
  #                      start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
  #                 control=ctrl, lower = c(0,0,0), upper = c(1,1,100))
  #
  # s = summary(model1)
  # coeff1 = s$parameters[,1]
  # r1 = sum((s$residuals)^2)

  ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)

  model1 <- nlxb(z ~ a1 / (1+exp(-a2*(tk-a3))),
                  start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
                 control = ctrl, lower = c(0,0,0), upper = c(1,1,100))

  coeff1 = model1$coefficients
  r1 = model1$ssquares


  # model2 <- nlsLM(z ~ a1 + exp(-a2*(tk+a3)),
  #                 start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
  #                 control=ctrl, lower = c(0,0,0), upper = c(1,1,100))
  # s2=summary(model2)
  # coeff2 = s2$parameters[,1]
  # r2 = sum((s2$residuals)^2)

  ctrl = list(phi=1, lamda = 0.0087, offset = 100, laminc=10, lamdec = 4)


  model2 <- nlxb(z ~ a1 + exp(-a2*(tk+a3)),
                  start=list(a1=guess[5], a2=guess[6], a3=guess[7]), data=df, trace=trace,
                 control = ctrl, lower = c(0,0,0), upper = c(1,1,100))
  coeff2 = model2$coefficients
  r2 = model2$ssquares


  # coeff2 = pracma::lsqcurvefit(myFun2, guess[5:7], x[!is.na(rate)], rate[!is.na(rate)])$x
  # r2 = pracma::lsqcurvefit(myFun2, guess[5:7], x[!is.na(rate)], rate[!is.na(rate)])$ssq
  #
  #
  # %                  figure;plot(x,rate,x,myFun1(coeff1,x),'r',x,myFun2(coeff2,x),'g--')

  # % myFun1 is more stable on a long term persepective
  # %               % If coeff2 have reached the upper boundaries, myFUn1 is
  # %               chosen
  if (r1<r2 || coeff2[1]>0.99 || coeff2[2]>4.9){
    lambdaGuess = coeff1
    lambdaFun = myFun1
  }else{
    lambdaGuess = coeff2
    lambdaFun = myFun2
  }


  guess[5:7] = lambdaGuess; #% update guess

  }
  return(list(guess=guess, lambdaFun = lambdaFun))

  # catch exceptionL
  # disp(exceptionL)
  # lambdaFun =  @(a,t) a(1)./(1+exp(-a(2)*(t-a(3))));
}
