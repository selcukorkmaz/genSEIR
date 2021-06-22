#' Plots for Epidemic Curves
#'
#' This function creates plots for reported and predicted active, recovered and death cases.

#' @param object a SEIQRDP result
#' @param active time histories of the quarantined/active cases
#' @param recovered time histories of the recovered cases
#' @param deaths Time histories of the deceased cases
#' @param title country name
#' @param params fitted parameters by fit_SEIQRDP function
#' @param checkRates if \code{TRUE} compares the fitted and calcualted death and recovered ratios through plots
#' @param ... other plot options
#'
#' @importFrom ggplot2  ggplot geom_line geom_point aes scale_colour_manual scale_fill_manual xlab ylab ggtitle labs
#'
#' @export plotSEIQRDP
#'
#' @return plots for epidemic curves: active cases, recovered and deaths
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @examples
#' \donttest{
#' start = "01/01/21"
#' finish = "04/01/21"
#' country = "Italy"
#' dt = 1
#' f=30
#'
#' covidData = getDataCOVID(start = start, finish = finish, country = country)
#' Recovered = covidData$tableRecovered
#' Deaths = covidData$tableDeaths
#' Confirmed = covidData$tableConfirmed
#'
#' if(nrow(Recovered) == 1){
#'   name = Recovered$CountryRegion
#' }else{
#'    name = paste0(Recovered$ProvinceState, " (",Recovered$CountryRegion,")")
#' }
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
#'
#'  res = SEIQRDP(alpha = params$alpha1, beta = params$beta1,
#'                gamma = params$gamma1, delta = params$delta1,
#'                lambda0 = params$Lambda1, kappa0 = params$Kappa1,
#'                Npop, E0, I0, Q0, R0, D0,lambdaFun = params$lambdaFun,
#'                kappaFun = params$kappaFun, tstart = start, tfinish = finish,
#'                dt = dt, f =f)
#'
#' p = plotSEIQRDP(res, Active, recovered, deaths, name, params, checkRates = TRUE)
#' print(p)
#' }
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}


plotSEIQRDP <- function(object, active, recovered, deaths, title = NULL, params, checkRates = FALSE, ...){

  realTime = as.Date(object$realTime)
  simTime = as.Date(object$simTime)

  pred = cbind.data.frame(S = object$susceptible, E = object$exposed, I = object$infectious,
             Q = object$quarantined, R = object$recovered, D = object$dead,
             P = object$insusceptible)

  actual = cbind.data.frame(active = as.numeric(active), recovered = as.numeric(recovered), deaths = as.numeric(deaths))
  color_fitted <- c("red", "black", "blue")
  color_filled = c("red", "black", "blue")

     fittedPlot <- ggplot(pred) +
     geom_line(aes(y = R, x=simTime, color = "Recovered")) +
     geom_line(aes(y = Q, x=simTime, color="Active")) +
     geom_line(aes(y = D, x=simTime, color="Deceased"))+
     # geom_line(aes(y = I, x=simTime, color="Infectious"))+
     # geom_line(aes(y = E, x=simTime, color="Exposed"))+
     scale_colour_manual(values=color_fitted) +
     xlab("Time") + ylab("Cases")+ggtitle(title)+
     geom_point(data = actual,  aes(y = recovered, x = realTime, fill = "Recovered"), size=2, color = "blue", shape=21, stroke=0)+
     geom_point(data = actual,  aes(y = active, x = realTime, fill = "Active"), size=2,color = "red", shape=21, stroke=0)+
     geom_point(data = actual,  aes(y = deaths, x = realTime, fill = "Deceased"), size=2,color="black",shape=21, stroke=0)+
     scale_fill_manual(values=color_filled)+
     labs(fill="Reported", color="Fitted")

     if(checkRates){

      checkRates = checkRates(object$realTime, active,recovered, deaths,
                              params$kappaFun, params$lambdaFun, params$Kappa1, params$Lambda1, object$dt)

     }

    return(list(fittedPlot, checkRates))

}

