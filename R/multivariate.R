#############################################
# AB lasso mv with random sample splitting
#############################################
#' AB-LASSO Estimator with Random Sample Splitting for Multivariate Models
#'
#' Implements the AB-LASSO estimation method for the multivariate model
#' \eqn{Y_{it} = \alpha_{i} + \gamma_{t} + \sum_{j=1}^{L} \beta_{j} Y_{i,t-j} + \theta_{0} D_{it} + \theta_{1} C_{i,t-1} + \varepsilon_{it}}, with random sample splitting. Note that \eqn{D_{it}} and \eqn{C_{it}} are predetermined with respect to \eqn{\varepsilon_{it}}.
#'
#' @param Y A `P` x `N` (number of time periods x number of individuals) matrix containing the outcome/response variable `Y`.
#' @param D A `P` x `N` (number of time periods x number of individuals) matrix containing the policy variable/treatment `D`.
#' @param C A list of `P` x `N` matrices containing other treatments and control variables.
#' @param lag The lag order of \eqn{Y_{it}} included in the covariates, default is `1`.
#' @param Kf The number of folds for K-fold cross-validation, with options being `2` or `5`, default is `2`.
#' @param nboot The number of random sample splits, default is `100`.
#' @param seed Seed for random number generation, default `202302`.
#'
#' @return A dataframe that includes the estimated coefficients (\eqn{\beta_{j}, \theta_{0}, \theta_{1}}), their standard errors, and T-statistics.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Use the Covid data
#' N = length(unique(covid_data$fips))
#' P = length(unique(covid_data$week))
#' Y = matrix(covid_data$logdc, nrow = P, ncol = N)
#' D = matrix(covid_data$dlogtests, nrow = P, ncol = N)
#' C = list()
#' C[[1]] = matrix(covid_data$school, nrow = P, ncol = N)
#' C[[2]] = matrix(covid_data$college, nrow = P, ncol = N)
#' C[[3]] = matrix(covid_data$pmask, nrow = P, ncol = N)
#' C[[4]] = matrix(covid_data$pshelter, nrow = P, ncol = N)
#' C[[5]] = matrix(covid_data$pgather50, nrow = P, ncol = N)
#'
#' results.kf2 <- ablasso_mv_ss(Y = Y, D = D, C = C, lag = 4, nboot = 2)
#' print(results.kf2)
#' results.kf5 <- ablasso_mv_ss(Y = Y, D = D, C = C, lag = 4, Kf = 5, nboot = 2)
#' print(results.kf5)
#' }
ablasso_mv_ss <- function(Y, D, C,
                          lag=1, Kf=2,
                          nboot=100, seed=202302) {
  Kf_list <- c(2, 5)
  if (!Kf %in% Kf_list) {
    stop("Invalid Kf type. Expected one of {2, 5}.")
  }

  # Check if Y and D are matrices or arrays of the same size
  if (!is.matrix(Y) || !is.matrix(D) || !all(dim(Y) == dim(D))) {
    stop("Y and D must be matrices or arrays of the same size.")
  }

  # Check each matrix in C to ensure it has the correct dimensions and is a matrix
  for (i in seq_along(C)) {
    if (!is.matrix(C[[i]])) {
      stop(sprintf("C[[%d]] is not a matrix.", i))
    }
    if (!all(dim(Y) == dim(C[[i]]))) {
      stop(sprintf("C[[%d]] must have the same dimensions as Y and D.", i))
    }
  }

  N = ncol(Y)
  P = nrow(Y)
  Y.t = diff(Y)
  D.t = diff(D)
  C.t = list()
  for(j in 1:length(C)){
    C.t[[j]] = diff(C[[j]])
  }

  W1 = list()
  for(j in 1:lag){
    W1[[j]] = Y.t[(lag-j+1):(nrow(Y.t)-j),]
  }
  W2 = D.t[(lag+1):nrow(D.t),]
  W4 = array(0, dim = c(P-lag-1,N,length(C)))
  for(j in 1:length(C)){
    W4[,,j] = C.t[[j]][lag:(nrow(C.t[[j]])-1),]
  }
  Q = array(0, dim = c(P,N,P-lag-1))
  for(t in (lag+2):P){
    Q[t,,t-lag-1] = rep(1, N)
  }
  Q.t = array(0, dim = c(P-1,N,P-lag-1))
  for(j in 1:(P-lag-1)){
    Q.t[,,j] = diff(Q[,,j])
  }
  W3 = array(0, dim = c(P-lag-1,N,P-lag-1))
  for(j in 1:(P-lag-1)){
    W3[,,j] = Q.t[(lag+1):nrow(Q.t[,,j]),,j]
  }

  ############### AB-LASSO-SS ################
  set.seed(seed)
  theta.hat.all = std.hat.all = numeric(0)
  vcv.all = list()
  for(ib in 1:nboot){
    foldid = rep.int(1:Kf, times=ceiling(N/Kf))[sample.int(N)] #fold IDs
    I = split(1:N, foldid)

    Z = list()
    for(t in (lag+2):P){
      Z[[t-lag-1]] = matrix(0, N, (1+lag+dim(W3)[3]+1+length(C)-1))
      y1 = y2 = rep(0, N)
      y3 = matrix(0, N, dim(W3)[3])
      y4 = matrix(0, N, dim(W4)[3])
      zz = matrix(0, N, (length(C)*(t-1)))
      x = matrix(0, N, ((t-2)+(t-1)+dim(zz)[2]))
      for(b in 1:length(I)){
        y1[-I[[b]]] = W1[[1]][t-lag-1,-I[[b]]]  #auxiliary sample
        y2[-I[[b]]] = W2[t-lag-1,-I[[b]]]
        y3[I[[b]],] = W3[t-lag-1,I[[b]],]       #main sample
        y4[I[[b]],] = W4[t-lag-1,I[[b]],]
        for(j in 1:length(C)){
          zz[I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(C[[j]][1:(t-1),I[[b]]])
          zz[-I[[b]],((j-1)*(t-1)+1):(j*(t-1))] = t(C[[j]][1:(t-1),-I[[b]]])
        }
        if(lag==1){if(t>lag+2){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}else{x[I[[b]],] = cbind(Y[1:(t-2),I[[b]]], t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}}
        if(lag>1){x[I[[b]],] = cbind(t(Y[1:(t-2),I[[b]]]), t(D[1:(t-1),I[[b]]]), zz[I[[b]],])}
        if(lag==1){if(t>lag+2){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}else{x[-I[[b]],] = cbind(Y[1:(t-2),-I[[b]]], t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}}
        if(lag>1){x[-I[[b]],] = cbind(t(Y[1:(t-2),-I[[b]]]), t(D[1:(t-1),-I[[b]]]), zz[-I[[b]],])}
        fit1_s2 = hdm::rlasso(x[-I[[b]],], y1[-I[[b]]])
        Z[[t-lag-1]][I[[b]],1] = stats::predict(fit1_s2, x[I[[b]],])
        j = 1
        while(j<lag){
          Z[[t-lag-1]][I[[b]],j+1] = W1[[j+1]][t-lag-1,I[[b]]]
          j = j + 1
        }
        fit2_s2 = hdm::rlasso(x[-I[[b]],], y2[-I[[b]]])
        Z[[t-lag-1]][I[[b]],lag+1] = stats::predict(fit2_s2, x[I[[b]],])
        Z[[t-lag-1]][I[[b]],(lag+1+1):(2+lag+dim(W4)[3]-1)] = y4[I[[b]],]
        Z[[t-lag-1]][I[[b]],(lag+2+dim(W4)[3]):(2+lag+dim(W4)[3]+dim(W3)[3]-1)] = y3[I[[b]],]
      }
    }

    theta.hat = matrix(0, (1+lag+dim(W3)[3]+1+length(C)-1), length(I))
    for(b in 1:length(I)){
      sum1 = matrix(0, P-lag-1+lag+1+1+length(C)-1, P-lag-1+lag+1+1+length(C)-1)
      sum2 = rep(0, P-lag-1+lag+1+1+length(C)-1)
      for(i in I[[b]]){
        for(t in (lag+2):P){
          j = 1
          XX = numeric(0)
          while(j<=lag){
            XX = c(XX, W1[[j]][t-lag-1,i])
            j = j + 1
          }
          sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,])) # lag of Y, D, C
          sum2 = sum2 + Z[[t-lag-1]][i,]*Y.t[t-1,i]
        }
      }
      sum1.inv = try(solve(sum1))
      theta.hat[,b] = sum1.inv%*%sum2
    }
    theta.hat = rowMeans(theta.hat)

    W4.all = numeric()
    for(j in 1:dim(W4)[3]){
      W4.all = cbind(W4.all, as.vector(W4[,,j]))
    }
    W3.all = numeric()
    for(j in 1:dim(W3)[3]){
      W3.all = cbind(W3.all, as.vector(W3[,,j]))
    }
    j = 1
    XX = numeric(0)
    while(j<=lag){
      XX = cbind(XX, as.vector(W1[[j]]))
      j = j + 1
    }
    vps.t = matrix(as.vector(Y.t[(lag+1):nrow(Y.t),]) - cbind(XX, as.vector(W2), W4.all, W3.all)%*%theta.hat, nrow = P-lag-1, ncol = N)
    mu = matrix(0, N, P-lag-1+lag+1+1+length(C)-1)
    for(i in 1:N){
      for(t in (lag+2):P){
        mu[i,] = mu[i,] + Z[[t-lag-1]][i,]*vps.t[t-lag-1,i]
      }
      mu[i,] = mu[i,]/(P-lag-1)
    }
    sum3 = sum4 = matrix(0, P-lag-1+lag+1+1+length(C)-1, P-lag-1+lag+1+1+length(C)-1)
    for(i in 1:N){
      for(t in (lag+2):P){
        sum3 = sum3 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])
      }
    }
    for(i in 1:N){
      for(t in (lag+2):(P-1)){
        sum4 = sum4 + (Z[[t-lag-1]][i,]*vps.t[t-lag-1,i] - mu[i,])%*%t(Z[[t-lag]][i,]*vps.t[t-lag,i] - mu[i,])
      }
    }
    Sigma = sum3 + sum4*(P-lag-2)/(P-lag-1) + t(sum4)*(P-lag-2)/(P-lag-1)
    sum1 = matrix(0, P-lag-1+lag+1+1+length(C)-1, P-lag-1+lag+1+1+length(C)-1)
    for(i in 1:N){
      for(t in (lag+2):P){
        j = 1
        XX = numeric(0)
        while(j<=lag){
          XX = c(XX, W1[[j]][t-lag-1,i])
          j = j + 1
        }
        sum1 = sum1 + Z[[t-lag-1]][i,]%*%t(c(XX, W2[t-lag-1,i], W4[t-lag-1,i,], W3[t-lag-1,i,]))
      }
    }
    sum1.inv = try(solve(sum1))
    std.hat = sqrt(diag(sum1.inv%*%Sigma%*%t(sum1.inv)))
    vcv = sum1.inv%*%Sigma%*%t(sum1.inv)

    theta.hat.all = rbind(theta.hat.all, theta.hat)
    std.hat.all = rbind(std.hat.all, std.hat)
    vcv.all[[ib]] = vcv

    message(paste(ib, "/", nboot))
  }

  num_coef <- lag + 1 + length(C)
  theta.hat.all.med <- matrixStats::colMedians(theta.hat.all)[1:num_coef] # the order of the coefficients: lags of Y, D_{t}, C_{t-1}, and time dummies
  std.hat.all.med <- matrixStats::colMedians(std.hat.all)[1:num_coef]  # std. errors
  stat.all.med <- matrixStats::colMedians(theta.hat.all)[1:num_coef]/matrixStats::colMedians(std.hat.all)[1:num_coef] # T-stat

  res <- rbind(estimates = theta.hat.all.med,
               std.hat = std.hat.all.med,
               stat = stat.all.med)

  results <- as.data.frame(res)

  # rename the columns
  for (i in 1:lag) {
    colnames(results)[i] <- paste("beta", i, sep = "_")
  }
  colnames(results)[lag+1] <- "theta_0"
  j = 1
  for (i in (lag+2):(num_coef)) {
    colnames(results)[i] <- paste("theta_1", j, sep = ".")
    j = j + 1
  }

  return(results)
}


