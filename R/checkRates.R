
checkRates <- function(time,Q,R,D,kappaFun,lambdaFun,kappa,lambda, dt=1){
    #' Check Rates
    #'
    #' This function compares the fitted and calcualted death and recovered ratios.
    #' The idea is to check whether the approximation of these ratios is appropriate

    #' @param time Time vector
    #' @param Q Time histories of the quarantined/active cases
    #' @param R Time histories of the recovered cases
    #' @param D Time histories of the deceased cases
    #' @param kappaFun Anonymous function approximating the death rate
    #' @param lambdaFun Anonymous function approximating the recovery rate
    #' @param kappa mortality rate
    #' @param lambda cure rate
    #' @param dt Time step, default is 1
    #'
    #' @export checkRates
    #'
    #' @return Creates plots for death rate and recovery rate. The idea is to check whether the approximation of these ratios is appropriate.
    #'
    #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
    #'
    #' @seealso \code{\link{SEIQRDP}} \code{\link{fit_SEIQRDP}}
    #'
    #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}
    #'


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
