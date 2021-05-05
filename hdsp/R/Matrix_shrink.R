# Shrinkage of the covariance matrix


#' For shrinkage of the covariance matrix (2014, JMVA, Strong convergence)
#'
#' It returns the estimator CovShrink(TM,SCM)=alfa*SCM+beta*TM,
#' where SCM should be the sample covariance matrix and TM is a shrinkage target.
#' In most cases, we take TM=I but could be any deterministic pos. def. matrix.
#' In order to use the estimator for portfolio weights one needs to invert it,
#' i.e., solve(CovShrink(TM,SCM)).
#' @param n the number of observations of x.
#' @param TM a shrinkage target.
#' @param SCM the sample covariance matrix.
#' @references \insertRef{BGP2014}{hdsp}
#' @examples
#' # Parameter setting
#' n<-5e2
#' c<-0.7
#' p<-c*n
#' mu <- rep(0, p)
#' Sigma <- RandCovMtrx(n=n, p=p, q=20.55)
#'
#' # Generating observations
#' X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))
#'
#' # Estimation
#' TM <- matrix(0, nrow=p, ncol=p)
#' diag(TM) <- 1
#' SCM <- Sigma_sample_estimator(X)
#' Sigma_shr <- CovShrinkBGP14(n=n, TM=TM, SCM=SCM)
#' Sigma_shr[1:6, 1:6]
#' @export
CovShrinkBGP14<-function(n, TM, SCM)
{
  a_1<-(1/n)*sum(diag(TM%*%TM))*(sum(diag(SCM)))^2
  a_2<-sum(diag(SCM%*%SCM))*sum(diag(TM%*%TM))-(sum(diag(SCM%*%TM)))^2
  alfa1<-(1-a_1/a_2)
  beta1<-sum(diag(SCM%*%TM))*(1-alfa1)/sum(diag(TM%*%TM))
  BSR<-alfa1*SCM+beta1*TM
  return(BSR)
}


#' Ledoit-Wolf Nonlinear shrinkage estimator (Annals of Statistics (2020))
#'
#' Returns the estimator of the covariance matrix.
#'
#' @param x a numeric matrix. Rows represent different variables, columns- observations.
#' @references \insertRef{LW2020}{hdsp}
#' @examples
#' n<-5e2
#' c<-0.7
#' p<-c*n
#' mu <- rep(0, p)
#' Sigma <- RandCovMtrx(n=n, p=p, q=20.55)
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


# Shrinkage of the inverse covariance matrix

#' For the shrinkage of the inverse covariance matrix
#'
#' iSCM=solve(S) with S sample covariance matrix and TM is again a target matrix,
#' for example TM=I. Thus, InvCovShrink(solve(S), TM) will return the estimator of
#' the inverse covariance matrix (no need to invert anymore). n- the number of observations of X.
#' @param TM the target matrix for shrinkage of the inverse covariance matrix
#' @param n the number of observations
#' @param p the number of variables (rows of the covariance matrix)
#' @param iSCM the inverse of the sample covariance matrix
#' @references \insertRef{BGP2016}{hdsp}
#' @export
InvCovShrinkBGP16<-function(n, p, TM, iSCM)
{
  a_1<-(1/n)*sum(diag(TM%*%TM))*(sum(diag(iSCM)))^2
  a_2<-sum(diag(iSCM%*%iSCM))*sum(diag(TM%*%TM))-(sum(diag(iSCM%*%TM)))^2
  alfa1<-(1-p/n-a_1/a_2)
  beta1<-sum(diag(iSCM%*%TM))*(1-p/n-alfa1)/sum(diag(TM%*%TM))
  BSR<-alfa1*iSCM+beta1*TM
  return(BSR)
}

