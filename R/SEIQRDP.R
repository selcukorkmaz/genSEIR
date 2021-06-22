#'Simulate generalized SEIR model
#'
#'This function simulates the time-histories of an epidemic outbreak using a generalized SEIR model
#'
#' @param alpha fitted protection rate
#' @param beta fitted  infection rate
#' @param gamma fitted  Inverse of the average latent time
#' @param delta fitted  rate at which people enter in quarantine
#' @param lambda0 fitted  cure rate
#' @param kappa0 fitted  mortality rate
#' @param Npop Total population of the sample
#' @param E0 Initial number of exposed cases
#' @param I0 Initial number of infectious cases
#' @param Q0 Initial number of quarantined cases
#' @param R0 Initial number of recovered cases
#' @param D0 Initial number of dead cases
#' @param lambdaFun anonymous function giving the time-dependant recovery rate
#' @param kappaFun anonymous function giving the time-dependant death rate
#' @param tstart start date
#' @param tfinish finish date
#' @param dt the time step. This oversamples time to ensure that the algorithm converges
#' @param f future predictions
#'
#' @export SEIQRDP
#'
#' @importFrom stats median
#'
#'@author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @return a list of predicted cases including susceptible, exposed,  infectious, quarantined, recovered, dead and insusceptible.
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
#'}
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @seealso \code{\link{fit_SEIQRDP}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

SEIQRDP <- function(alpha, beta, gamma, delta, lambda0, kappa0,
                    Npop, E0, I0, Q0, R0, D0, lambdaFun, kappaFun,
                    tstart, tfinish, dt = 1/24, f=0){

realTime = seq(as.Date(tstart, format = "%m/%d/%y"), as.Date(tfinish, format = "%m/%d/%y"), by = "1 day")
simTime = seq(as.Date(realTime[1]),as.Date(realTime[length(realTime)])+f, dt)
N = length(simTime)
t = (1:N-1)*dt

##Initial conditions
Y = matrix(0,7,N)
Y[1,1] = as.numeric(Npop-Q0-E0-R0-D0-I0)
Y[2,1] = as.numeric(E0)
Y[3,1] = as.numeric(I0)
Y[4,1] = as.numeric(Q0)
Y[5,1] = as.numeric(R0)
Y[6,1] = as.numeric(D0)

if (round(sum(Y[,1])-Npop)!=0){
  stop('the sum must be zero because the total population (including the deads) is assumed constant')
}

# Computes the seven states
dt = median(diff(t))
lambda1 = lambdaFun(lambda0,t)
kappa1 = kappaFun(kappa0,t)

# ODE resolution
for(ii in 1:(N-1)){
    A = getA(alpha,gamma,delta,lambda1[ii],kappa1[ii])
    SI = Y[1,ii]*Y[3,ii]
    K = matrix(0,7)
    K[1:2,1] = rbind(-beta/Npop,beta/Npop)*SI
    Y[,ii+1] = RK4(Y[,ii],A,K,dt)
}

 # Y = round(Y)
 # Write the outputs
  S = Y[1,1:N]
  E = Y[2,1:N]
  I = Y[3,1:N]
  Q = Y[4,1:N]
  R = Y[5,1:N]
  D = Y[6,1:N]
  P = Y[7,1:N]



  return(list(susceptible = S,exposed = E, infectious = I,quarantined = Q,
              recovered = R, dead = D, insusceptible = P, realTime = realTime, simTime = simTime, dt = dt))

}





