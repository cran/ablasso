#############################################
# data generating process
#############################################
#' Generate a Dataset for Simulations
#'
#' Generates data according to the following process:
#' \eqn{Y_{it} = \alpha_{i} + \gamma_{t} + \theta_{1} Y_{i,t-1} + \theta_{2} D_{it} + \varepsilon_{it}} and
#' \eqn{D_{it} = \rho D_{i,t-1} + v_{i,t}}.
#' Note that \eqn{D_{it}} is predetermined with respect to \eqn{\varepsilon_{it}}.
#'
#' @param N An integer specifying the number of individuals.
#' @param P An integer specifying the number of time periods.
#' @param sigma_alpha Standard deviation for the normal distribution from which the individual effect `alpha` is drawn; default is 1.
#' @param sigma_gamma Standard deviation for the normal distribution from which the time effect `gamma` is drawn; default is 1.
#' @param sigma_eps.d Standard deviation for the error term associated with the policy variable/treatment (`D`); default is `1`.
#' @param sigma_eps.y Standard deviation for the error term associated with the outcome/response variable (`Y`); default is `1`.
#' @param cov_eps Covariance between error terms of `Y` and `D`, default `0.5`.
#' @param rho Autocorrelation coefficient for `D` across time, default `0.5`.
#' @param theta Regression Coefficients for univariate AR(1) dynamic panal, default `c(0.8, 1)`.
#' @param seed Seed for random number generation, default `202304`.
#'
#' @return A list of two `P` x `N` matrices named `Y` (outcome/response variable) and `D` (policy variable/treatment).
#'
#' @export
#'
#' @examples
#' # Generate data using default parameters
#' data1 <- generate_data(N = 300, P = 40)
#' str(data1)
#'
#' data2 <- generate_data(N = 500, P = 20)
#' str(data2)
generate_data <- function(N, P,
                          sigma_alpha=1, sigma_gamma=1,
                          sigma_eps.d=1, sigma_eps.y=1,
                          cov_eps=0.5, rho=0.5, theta=c(0.8,1),
                          seed=202304) {

  # Initialize a list to store results
  data = list()
  set.seed(seed)

  # Generate alpha, gamma, eps, ...
  alpha = matrix(0, nrow = 2, ncol = N)
  gamma = matrix(0, nrow = 2, ncol = P+10)
  eps.y = eps.d = array(0, dim = c(2, P+10, N))
  eps = array(0, dim = c(2, P+10, N, 2))

  for(r in 1:2){
    alpha[r,] = stats::rnorm(N, 0, sigma_alpha)
    gamma[r,] = stats::rnorm(P+10, 0, sigma_gamma)
    eps.y[r,,] = mvtnorm::rmvnorm(P+10, rep(0,N), diag(N)*sigma_eps.y)
    eps.d[r,,] = mvtnorm::rmvnorm(P+10, rep(0,N), diag(N)*sigma_eps.d)
    for(i in 1:N){
      eps[r,,i,] = mvtnorm::rmvnorm(P+10, rep(0,2), cbind(c(sigma_eps.y, cov_eps), c(cov_eps, sigma_eps.d)))
    }
  }

  alpha = alpha[1,]
  gamma = gamma[1,]
  eps.y = eps.y[1,,]
  eps.d = eps.d[1,,]
  eps = eps[1,,,]
  for(i in 1:N){
    eps.y[,i] = eps[,i,1]
    eps.d[2:(P+10),i] = eps[1:(P+10-1),i,2]
  }

  # Initialize Y and D
  Y = D = matrix(0, nrow = P+10, ncol = N)
  for(t in 2:(P+10)){
    D[t,] = rho*D[t-1,] + eps.d[t,]
    Y[t,] = theta[1]*Y[t-1,] + theta[2]*D[t,] + alpha + eps.y[t,] + gamma[t]
  }
  Y = Y[-(1:10),]
  D = D[-(1:10),]

  # Return a list that contains Y and D
  data = list(Y = Y, D = D)
  return(data)
}


#############################################
# check data_structure
#############################################
check_data_structure <- function(Y, D) {

  # Check if Y and D are matrices or arrays of the same size
  if (!is.matrix(Y) || !is.matrix(D) || !all(dim(Y) == dim(D))) {
    stop("Y and D must be matrices or arrays of the same size.")
  }

}


#############################################
# AB-LASSO without sample splitting
#############################################
#' AB-LASSO Estimator Without Sample Splitting
#'
#' Implements the AB-LASSO estimation method for the univariate model \eqn{Y_{it} = \alpha_{i} + \gamma_{t} + \theta_{1} Y_{i,t-1} + \theta_{2} D_{it} + \varepsilon_{it}}, without sample splitting. Note that \eqn{D_{it}} is predetermined with respect to \eqn{\varepsilon_{it}}.
#'
#' @param Y A `P` x `N` (number of time periods x number of individuals) matrix containing the outcome/response variable `Y`.
#' @param D A `P` x `N` (number of time periods x number of individuals) matrix containing the policy variable/treatment `D`.
#'
#' @return A list with three elements:
#' \itemize{
#' \item theta.hat: Estimated coefficients.
#' \item std.hat: Estimated Standard errors.
#' \item stat: T-Statistics.
#' }
#' @export
#'
#' @examples
#' # Generate data
#' data1 <- generate_data(N = 300, P = 40)
#'
#' # You can use your own data by providing matrices `Y` and `D`
#' results <- ablasso_uv(Y = data1$Y, D = data1$D)
#' print(results)
ablasso_uv <- function(Y, D) {

  check_data_structure(Y, D)

  results = numeric(0)
  output = list()
  theta.hat = matrix(0, 1, 2)

  # Identify N, P and rep for each list
  N = dim(Y)[2]
  P = dim(Y)[1]

  Y.t = diff(Y) - rowMeans(diff(Y))
  D.t = diff(D) - rowMeans(diff(D))
  W1 = Y.t[1:(nrow(Y.t)-1),]
  W2 = D.t[2:nrow(D.t),]
  Z = list()
  for(t in 3:P){
    y1 = W1[t-2,]
    y2 = W2[t-2,]
    if(t>3){x = cbind(t(Y[1:(t-2),]), t(D[1:(t-1),]))}else{x = cbind(Y[1:(t-2),], t(D[1:(t-1),]))}
    fit1 = hdm::rlasso(x, y1)
    fit2 = hdm::rlasso(x, y2)
    Z[[t-2]] = cbind(stats::predict(fit1), stats::predict(fit2))
  }
  sum1 = matrix(0, 2, 2)
  sum2 = rep(0, 2)
  for(i in 1:N){
    for(t in 3:P){
      sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
      sum2 = sum2 + Z[[t-2]][i,]*Y.t[t-1,i]
    }
  }
  # calculate theta.hat
  theta.hat = solve(sum1)%*%sum2
  vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1),as.vector(W2))%*%theta.hat, nrow = P-2, ncol = N)
  mu = matrix(0, N, 2)
  for(i in 1:N){
    for(t in 3:P){
      mu[i,] = mu[i,] + Z[[t-2]][i,]*vps.t[t-2,i]
    }
    mu[i,] = mu[i,]/(P-2)
  }
  sum3 = sum4 = matrix(0, 2, 2)
  for(i in 1:N){
    for(t in 3:P){
      sum3 = sum3 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])
    }
  }
  for(i in 1:N){
    for(t in 3:(P-1)){
      sum4 = sum4 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
    }
  }
  Sigma = sum3 + sum4*(P-3)/(P-2) + t(sum4)*(P-3)/(P-2)

  std.hat = sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))*N*(P-2)))
  stat = sqrt(N*(P-2))*theta.hat/std.hat

  output = list(theta.hat = theta.hat,
                std.hat = std.hat,
                stat = stat)
  return(output)
}


#############################################
# random sample splitting
#############################################
#' AB-LASSO Estimator with Random Sample Splitting
#'
#' Implements the AB-LASSO estimation method for the univariate model \eqn{Y_{it} = \alpha_{i} + \gamma_{t} + \theta_{1} Y_{i,t-1} + \theta_{2} D_{it} + \varepsilon_{it}}, incorporating random sample splitting. Note that \eqn{D_{it}} is predetermined with respect to \eqn{\varepsilon_{it}}.
#'
#' @param Y A `P` x `N` (number of time periods x number of individuals) matrix containing the outcome/response variable variable `Y`.
#' @param D A `P` x `N` (number of time periods x number of individuals) matrix containing the policy variable/treatment `D`.
#' @param nboot The number of random sample splits, default is `100`.
#' @param Kf The number of folds for K-fold cross-validation, with options being `2` or `5`, default is `2`.
#' @param seed Seed for random number generation, default `202304`.
#'
#' @return A list with three elements:
#' \itemize{
#' \item theta.hat: Estimated coefficients.
#' \item std.hat: Estimated Standard errors.
#' \item stat: T-Statistics.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Generate data
#' data1 <- generate_data(N = 300, P = 40)
#'
#' # You can use your own data by providing matrices `Y` and `D`
#' results.ss <- ablasso_uv_ss(Y = data1$Y, D = data1$D, nboot = 2)
#' print(results.ss)
#'
#' results.ss2 <- ablasso_uv_ss(Y = data1$Y, D = data1$D, nboot = 2, Kf = 5)
#' print(results.ss2)
#' }
ablasso_uv_ss <- function(Y, D,
                       nboot=100,
                       Kf=2,
                       seed=202304) {

  check_data_structure(Y, D)
  Kf_list <- c(2, 5)
  if (!Kf %in% Kf_list) {
    stop("Invalid Kf type. Expected one of {2, 5}.")
  }

  set.seed(seed)
  output = list()

  N = dim(Y)[2]
  P = dim(Y)[1]

  theta.hat.ss = std.hat.ss = numeric(0)
  for(ib in 1:nboot){
    foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs
    I = split(1:N, foldid)

    #### Kf == 5 ####
    if(Kf==5){
      W1.all = W2.all = Y.t.all = list()
      for(b in 1:length(I)){
        Y.t = diff(Y)
        D.t = diff(D)
        Y.t[,I[[b]]] = Y.t[,I[[b]]] - rowMeans(Y.t[,I[[b]]])
        D.t[,I[[b]]] = D.t[,I[[b]]] - rowMeans(D.t[,I[[b]]])
        Y.t[,-I[[b]]] = Y.t[,-I[[b]]] - rowMeans(Y.t[,-I[[b]]])
        D.t[,-I[[b]]] = D.t[,-I[[b]]] - rowMeans(D.t[,-I[[b]]])
        W1 = Y.t[1:(nrow(Y.t)-1),]
        W2 = D.t[2:nrow(D.t),]
        W1.all[[b]] = W1
        W2.all[[b]] = W2
        Y.t.all[[b]] = Y.t
      }
      #### Kf == 2 ####
    }else{
      Y.t = diff(Y)
      D.t = diff(D)
      for(b in 1:length(I)){
        Y.t[,I[[b]]] = Y.t[,I[[b]]] - rowMeans(Y.t[,I[[b]]])
        D.t[,I[[b]]] = D.t[,I[[b]]] - rowMeans(D.t[,I[[b]]])
      }
      W1 = Y.t[1:(nrow(Y.t)-1),]
      W2 = D.t[2:nrow(D.t),]
    }

    Z = list()
    for(t in 3:P){
      Z[[t-2]] = matrix(0, N, 2)
      y1 = y2 = rep(0, N)
      x = matrix(0, N, ((t-2)+(t-1)))
      for(b in 1:length(I)){
        #### Kf == 5 ####
        if(Kf==5){
          W1 = W1.all[[b]]
          W2 = W2.all[[b]]
          Y.t = Y.t.all[[b]]
        }
        y1[-I[[b]]] = W1[t-2,-I[[b]]]  #s2 - auxiliary sample
        y2[-I[[b]]] = W2[t-2,-I[[b]]]
        if(t>3){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]))}else{x[I[[b]],] = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]))}
        if(t>3){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]))}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]))}
        fit1_s2 = hdm::rlasso(x[-I[[b]],], y1[-I[[b]]])
        fit2_s2 = hdm::rlasso(x[-I[[b]],], y2[-I[[b]]])
        Z[[t-2]][I[[b]],] = cbind(stats::predict(fit1_s2, x[I[[b]],]), stats::predict(fit2_s2, x[I[[b]],]))
      }
    }
    theta.est = matrix(0, 2, length(I))
    for(b in 1:length(I)){
      #### Kf == 5 ####
      if(Kf==5){
        W1 = W1.all[[b]]
        W2 = W2.all[[b]]
        Y.t = Y.t.all[[b]]
      }
      sum1 = matrix(0, 2, 2)
      sum2 = rep(0, 2)
      for(i in I[[b]]){
        for(t in 3:P){
          sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
          sum2 = sum2 + Z[[t-2]][i,]*Y.t[t-1,i]
        }
      }
      theta.est[,b] = solve(sum1)%*%sum2
    }
    theta.hat.ss = rbind(theta.hat.ss, rowMeans(theta.est))
    vps.t = matrix(as.vector(Y.t[2:nrow(Y.t),]) - cbind(as.vector(W1[,]),as.vector(W2[,]))%*%rowMeans(theta.est), nrow = P-2, ncol = N)
    mu = matrix(0, N, 2)
    for(i in 1:N){
      for(t in 3:P){
        mu[i,] = mu[i,] + Z[[t-2]][i,]*vps.t[t-2,i]
      }
      mu[i,] = mu[i,]/(P-2)
    }
    sum3 = sum4 = matrix(0, 2, 2)
    for(i in 1:N){
      for(t in 3:P){
        sum3 = sum3 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])
      }
    }
    for(i in 1:N){
      for(t in 3:(P-1)){
        sum4 = sum4 + (Z[[t-2]][i,]*vps.t[t-2,i] - mu[i,])%*%t(Z[[t-1]][i,]*vps.t[t-1,i] - mu[i,])
      }
    }
    Sigma = sum3 + sum4*(P-3)/(P-2) + t(sum4)*(P-3)/(P-2)
    sum1 = matrix(0, 2, 2)
    for(i in 1:N){
      for(t in 3:P){
        sum1 = sum1 + Z[[t-2]][i,]%*%t(c(W1[t-2,i],W2[t-2,i]))
      }
    }
    std.hat.ss = rbind(std.hat.ss, sqrt(diag(solve(sum1)%*%Sigma%*%t(solve(sum1))*N*(P-2))))
  }
  theta.hat = matrixStats::colMedians(theta.hat.ss)
  std.hat = matrixStats::colMedians(std.hat.ss)
  sd = matrixStats::colSds(theta.hat.ss)
  stat = sqrt(N*(P-2))*theta.hat/std.hat

  output = list(theta.hat = theta.hat,
                std.hat = std.hat,
                stat = stat)
  return(output)
}


