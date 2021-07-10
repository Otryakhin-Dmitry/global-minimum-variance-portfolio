# Shrinkage of the covariance matrix


#' Linear shrinkage estimator of the covariance matrix of \insertCite{BGP2014}{HDShOP}
#'
#' The optimal linear shrinkage estimator of the covariance matrix that minimizes the Frobenius norm:
#' \deqn{\hat{\Sigma}_{OLSE} = \hat{\alpha} \hat{\Sigma} + \hat{\beta} \Sigma_0,}
#' where \eqn{\hat{\alpha}} and \eqn{\hat{\beta}} are optimal shrinkage intensities
#' given in Eq. 4.3 and 4.4 of \insertCite{BGP2014}{HDShOP}. \eqn{\hat{\Sigma}}
#' is the sample covariance matrix (SCM) and \eqn{\Sigma_0} is a positive definite
#' symmetric matrix used as the target (TM), for example, \eqn{\frac{1}{p} I}.
#'
#' @param n sample size.
#' @param TM the target matrix for the shrinkage estimator.
#' @param SCM sample covariance matrix.
#'
#' @return a list containing an object of class matrix (S) and the estimated shrinkage
#' intensities \eqn{\alpha} and \eqn{\beta}.
#' @references \insertAllCited{}
#' @examples
#' # Parameter setting
#' n<-5e2
#' c<-0.7
#' p<-c*n
#' mu <- rep(0, p)
#' Sigma <- RandCovMtrx(p=p)
#'
#' # Generating observations
#' X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))
#'
#' # Estimation
#' TM <- matrix(0, nrow=p, ncol=p)
#' diag(TM) <- 1
#' SCM <- Sigma_sample_estimator(X)
#' Sigma_shr <- CovShrinkBGP14(n=n, TM=TM, SCM=SCM)
#' Sigma_shr$S[1:6, 1:6]
#' @export
CovShrinkBGP14<-function(n, TM, SCM)
{
  a_1<-(1/n)*sum(diag(TM%*%TM))*(sum(diag(SCM)))^2
  a_2<-sum(diag(SCM%*%SCM))*sum(diag(TM%*%TM))-(sum(diag(SCM%*%TM)))^2
  alfa1<-(1-a_1/a_2)
  beta1<-sum(diag(SCM%*%TM))*(1-alfa1)/sum(diag(TM%*%TM))
  BSR<-alfa1*SCM+beta1*TM
  list(S=BSR, alpha=alfa1, beta=beta1)
}


#' nonlinear shrinkage estimator of the covariance matrix  of Ledoit  and Wolf (2020, AoS)
#'
#' The nonlinear shrinkage estimator of the covariance matrix, that minimizes the
#' minimum variance loss functions as defined in Eq 2.1 of Ledoit and Wolf (2020, AoS).
#'
#' @inheritParams CovarEstim
#'
#' @return an object of class matrix
#' @references \insertRef{LW2020}{HDShOP}
#' @examples
#' n<-5e2
#' c<-0.7
#' p<-c*n
#' mu <- rep(0, p)
#' Sigma <- RandCovMtrx(p=p)
#'
#' X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))
#' Sigma_shr <- nonlin_shrinkLW(X)
#' @export
nonlin_shrinkLW = function(x){
  # the original version suggested that p is # of columns
  p = nrow(x)
  n = ncol(x)
  sampleC = Sigma_sample_estimator(x)
  eig = eigen(sampleC)
  u = eig$vectors[,p:1]
  lambda = rev(eig$values)
  lambda = lambda[max(1, p-n+1):p]
  L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
  h = n^(-1/3)
  H = h * t(L)
  x <- (L - t(L)) / H # This is a different x than before
  ftilde = (3/4/sqrt(5)) * rowMeans(pmax(1-x^2/5, 0) / H)
  Hftemp = (-3/10/pi) * x + (3/4/sqrt(5)/pi) * (1 - x^2./5) * log(abs((sqrt(5) - x)/(sqrt(5) + x)))
  Hftemp[abs(x) == sqrt(5)] = (-3/10/pi) * x[abs(x) == sqrt(5)]
  Hftilde = rowMeans(Hftemp / H)

  if(p<=n){
    dtilde = lambda / ((pi*(p/n)*lambda*ftilde)^2 + (1-(p/n)-pi*(p/n)*lambda*Hftilde)^2);
  }else{
    Hftilde0 = (1/pi)*(3/10/h^2+3/4/sqrt(5)/h*(1-1/5/h^2) *log((1+sqrt(5)*h)/(1-sqrt(5)*h))) * mean(1/lambda)
    dtilde0=1/(pi*(p-n)/n*Hftilde0)
    dtilde1 = lambda/(pi^2*lambda^2 * (ftilde^2+Hftilde^2))
    dtilde = c(dtilde0 * rep(1, p-n), dtilde1)
  }
  u %*% diag(dtilde) %*% t(u)
} # analytical nonlinear shrinkage


#' Linear shrinkage estimator of the inverse covariance matrix  of \insertCite{BGP2016}{HDShOP}
#'
#' The optimal linear shrinkage estimator of the inverse covariance (precision)
#' matrix that minimizes the Frobenius norm is given by:
#' \deqn{\hat{\Pi}_{OLSE} = \hat{\alpha} \hat{\Pi} + \hat{\beta} \Pi_0,}
#' where \eqn{\hat{\alpha}} and \eqn{\hat{\beta}} are optimal shrinkage intensities
#' given in Eq. 4.4 and 4.5 of \insertCite{BGP2016}{HDShOP}. \eqn{\hat{\Pi}} is
#' the sample inverse covariance matrix (iSCM) and \eqn{\Pi_0} is a positive definite
#' symmetric matrix used as the target (TM), for example, I.
#'
#' @param TM the target matrix for the shrinkage estimator
#' @param n the number of observations
#' @param p the number of variables (rows of the covariance matrix)
#' @param iSCM the inverse of the sample covariance matrix
#'
#' @return a list containing an object of class matrix (S) and the estimated shrinkage
#' intensities \eqn{\alpha} and \eqn{\beta}.
#' @references \insertAllCited{}
#' @examples
#' # Parameter setting
#' n<-5e2
#' c<-0.7
#' p<-c*n
#' mu <- rep(0, p)
#' Sigma <- RandCovMtrx(p=p)
#'
#' # Generating observations
#' X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))
#'
#' # Estimation
#' TM <- matrix(0, nrow=p, ncol=p)
#' diag(TM) <- 1
#' iSCM <- solve(Sigma_sample_estimator(X))
#' Sigma_shr <- InvCovShrinkBGP16(n=n, p=p, TM=TM, iSCM=iSCM)
#' Sigma_shr$S[1:6, 1:6]
#' @export
InvCovShrinkBGP16<-function(n, p, TM, iSCM)
{
  a_1<-(1/n)*sum(diag(TM%*%TM))*(sum(diag(iSCM)))^2
  a_2<-sum(diag(iSCM%*%iSCM))*sum(diag(TM%*%TM))-(sum(diag(iSCM%*%TM)))^2
  alfa1<-(1-p/n-a_1/a_2)
  beta1<-sum(diag(iSCM%*%TM))*(1-p/n-alfa1)/sum(diag(TM%*%TM))
  BSR<-alfa1*iSCM+beta1*TM
  list(S=BSR, alpha=alfa1, beta=beta1)
}

