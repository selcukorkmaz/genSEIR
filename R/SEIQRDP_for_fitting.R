modelFun = function(Y,A,K) {A %*% Y + K}

SEIQRDP_for_fitting <- function(par, t, t0){


  #'Nested Function
  #'
  #' @param par fitted protection rate
  #' @param t fitted  infection rate
  #' @param t0 fitted  inverse of the average latent time
  #'
  #' @importFrom pracma interp1

  #' @author Selçuk Korkmaz, \email{selcukorkmaz@gmail.com}
  #'
  #' @seealso \code{\link{fit_SEIQRDP}}
  #'
  #' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

    # % I simply rename the inputs
    alpha = abs(par[1])
    beta = abs(par[2])
    gamma = abs(par[3])
    delta = abs(par[4])
    lambda0 = abs(par[5:7])
    kappa0 = abs(par[8:10])


    # %% Initial conditions
    N = length(t)
    Y = data.frame(matrix(0,7,N)) #  There are seven different states
    Y[2,1] = E0;
    Y[3,1] = I0;
    Y[4,1] = Q[1];

    if (!is.null(R)){
    Y[5,1] = R[1]
    Y[1,1] = Npop-Q[1]-R[1]-D[1]-E0-I0
    }else{
      Y[1,1] = Npop-Q[1]-D[1]-E0-I0
    }
    Y[6,1] = D[1]

    if (round(sum(Y[,1])-Npop)!=0){
      stop('the sum must be zero because the total population including the deads) is assumed constant');
    }
    # %%

    kappa = kappaFun(kappa0,t)
    lambda = lambdaFun(lambda0,t)
    # % Very large recovery rate should not occur but can lead to
    # % numerical errors.

    if (length(lambda[lambda > 10])>0) {warning('lambda is abnormally high')}

    # % (ODE resolution
    for(ii in 1:(N-1)){
      A = getA(alpha,gamma,delta,lambda[ii],kappa[ii])
      SI = Y[1,ii]*Y[3,ii]
      F = matrix(0,7)
      F[1:2,1] = rbind(-beta/Npop,beta/Npop) %*% SI
      Y[,ii+1] = RK4(Y[,ii],A,F,dt)
    }

    Q1 = Y[4,1:N]
    R1 = Y[5,1:N]
    D1 = Y[6,1:N]

    Q1 = interp1(t,as.numeric(Q1[1,]),t0);
    R1 = interp1(t,as.numeric(R1[1,]),t0);
    D1 = interp1(t,as.numeric(D1[1,]),t0);

    if (!is.null(R)){
      output = cbind(Q1,R1,D1)
     colnames(output) = c("Q1", "R1", "D1")
    }else{
      output = cbind(Q1+R1, D1)
      colnames(output) = c("Q1", "D1")

    }

   return(output)
}
