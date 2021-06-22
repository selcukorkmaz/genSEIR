#' Check Rates
#'
#' This function compares the fitted and calculated death and recovered ratios.
#' The idea is to check whether the approximation of these ratios is appropriate.

#' @param time time vector
#' @param Q time histories of the quarantined/active cases
#' @param R time histories of the recovered cases
#' @param D time histories of the deceased cases
#' @param kappaFun anonymous function approximating the death rate
#' @param lambdaFun anonymous function approximating the recovery rate
#' @param kappa mortality rate
#' @param lambda cure rate
#' @param dt a time step, default is 1/24. This oversample time to ensure that the algorithm converges.
#'
#' @export checkRates
#'
#' @return plots for death rate and recovery rate
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

checkRates <- function(time, Q, R, D, kappaFun, lambdaFun, kappa,
                       lambda, dt = 1/24){

    Q = as.numeric(Q)
    R = as.numeric(R)
    D = as.numeric(D)

    rateD = (diff(D)/as.numeric(diff((time-time[1]))))/Q[2:length(Q)]
    rateD[abs(rateD)>3] = NA
    rateD[abs(rateD)<0] = NA

    if (!is.null(R)){
       rateR = (diff(R)/as.numeric(diff((time-time[1]))))/Q[2:length(Q)]
      rateR[abs(rateR)>3] = NA
    }

    x = as.numeric(time[2:length(time)]-time[1])
    x1 = seq(x[1], x[length(x)], dt)

    if (!is.null(R)){

        plot_death_rate <- ggplot() +
        geom_line(aes(y = kappaFun(kappa,x1), x=x1, color = "Fitted"))+
        scale_colour_manual(values="red") +
        geom_point(aes(y = rateD, x = x, fill = "Measured"),size=2, shape=21, color ="black", stroke=0)+
        scale_fill_manual(values="black")+
        labs(fill="", color="")+
        xlab("Time (days)") + ylab("Death rate")


        plot_recovery_rate <- ggplot() +
        geom_line(aes(y = lambdaFun(lambda,x1), x=x1, color = "Fitted"))+
        scale_colour_manual(values="blue") +
        geom_point(aes(y = rateR, x = x, fill = "Measured"),size=2, shape=21, color ="black", stroke=0)+
        scale_fill_manual(values="black")+
        labs(fill="", color="")+
        xlab("Time (days)") + ylab("Recovery rate")

        return(list(plot_death_rate, plot_recovery_rate))


    }else{

      plot_death_rate <- ggplot() +
        geom_line(aes(y = kappaFun(kappa,x1), x=x1, color = "Fitted"))+
        scale_colour_manual(values="red") +
        geom_point(aes(y = rateD, x = x, fill = "Measured"),size=2, shape=21, color ="black", stroke=0)+
        scale_fill_manual(values="black")+
        labs(fill="", color="")+
        xlab("Time (days)") + ylab("Pseudo-death rate")

      return(list(plot_death_rate))

    }
}
