SEIQRDP <- function(alpha,beta,gamma,delta,lambda0,kappa0,
                    Npop,E0,I0,Q0,R0,D0,lambdaFun,kappaFun,
                    tstart, tfinish, dt = 1, f=0){


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
  #' @param dt time step
  #' @param f future predictions
  #'
  #' @export SEIQRDP


  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}


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

  list=NA
}





