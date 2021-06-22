#'Fitted Results for SEIQRDP
#'
#' @param par initial guess parameters
#' @param t historical time vector
#' @param t0 target time vector
#'
#' @importFrom pracma interp1
#'
#' @return a data frame for fitted quarantined, recovered and deaths
#'
#' @author Selcuk Korkmaz, \email{selcukorkmaz@gmail.com}
#'
#' @seealso \code{\link{fit_SEIQRDP}}  \code{\link{RK4}}
#'
#' @references Peng, L., Yang, W., Zhang, D., Zhuge, C., Hong, L. 2020. “Epidemic analysis of COVID-19 in China by dynamical modeling”, arXiv preprint arXiv:2002.06563.
#' @references \url{https://www.mathworks.com/matlabcentral/fileexchange/74545-generalized-seir-epidemic-model-fitting-and-computation}

SEIQRDP_for_fitting <- function(par, t, t0){

    alpha = abs(par[1])
    beta = abs(par[2])
    gamma = abs(par[3])
    delta = abs(par[4])
    lambda0 = abs(par[5:7])
    kappa0 = abs(par[8:10])

    N = length(t)
    Y = data.frame(matrix(0,7,N))
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

    kappa = kappaFun(kappa0,t)
    lambda = lambdaFun(lambda0,t)

    if (length(lambda[lambda > 10])>0) {warning('lambda is abnormally high')}

    for(ii in 1:(N-1)){
      A = getA(alpha, gamma, delta, lambda[ii], kappa[ii])
      SI = Y[1,ii]*Y[3,ii]
      F = matrix(0,7)
      F[1:2,1] = rbind(-beta/Npop, beta/Npop) %*% SI
      Y[,ii+1] = RK4(Y[,ii], A, F, dt)
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

