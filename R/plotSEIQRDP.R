

plotSEIQRDP <- function(object, active, recovered, deaths, title = NULL, params, checkRates = FALSE, ...){

  #' Plot active, recovered and deaths
  #'
  #' This function computes the matrix A that is found in: dY/dt = A*Y + F

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
  #' @export plotSEIQRDP

  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}


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

