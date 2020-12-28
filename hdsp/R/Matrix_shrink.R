# Shrinkage of the covariance matrix


#' For shrinkage of the covariance matrix (2014, JMVA, Strong convergence)
#'
#' It returns the estimator CovShrink(TM,SCM)=alfa*SCM+beta*TM,
#' where SCM should be the sample covariance matrix and TM is a shrinkage target.
#' In most cases, we took TM=I but could be any deterministic pos. def. matrix.
#' In order to use the estimator for portfolio weights one needs to invert it,
#' i.e., solve(CovShrink(TM,SCM)).
CovShrinkBGP<-function(TM, SCM)
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
nonlin_shrinkLW = function(X){
  # the original version suggested that p is # of columns
  p = nrow(X)
  n = ncol(X)
  sampleC = Sigma_sample_estimator(X)
  eig = eigen(sampleC)
  u = eig$vectors[,p:1]
  lambda = rev(eig$values)
  lambda = lambda[max(1, p-n+1):p]
  L = matrix(rep(lambda, min(p, n)), nrow = length(lambda))
  h = n^(-1/3)
  H = h * t(L)
  x = (L - t(L)) / H
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

#' For the shrinkage of the inverse covariance matrix (2016, JMVA, Direct shrinkage)
#'
#' iSCM=solve(S) with S sample covariance matrix and TM is again a target matrix,
#' for example TM=I. Thus, InvCovShrink(solve(S), TM) will return the estimator of
#' the inverse covariance matrix (no need to invert anymore).
InvCovShrinkBGP<-function(TM,iSCM)
{
  a_1<-(1/n)*sum(diag(TM%*%TM))*(sum(diag(iSCM)))^2
  a_2<-sum(diag(iSCM%*%iSCM))*sum(diag(TM%*%TM))-(sum(diag(iSCM%*%TM)))^2
  alfa1<-(1-p/n-a_1/a_2)
  beta1<-sum(diag(iSCM%*%TM))*(1-p/n-alfa1)/sum(diag(TM%*%TM))
  BSR<-alfa1*iSCM+beta1*TM
  return(BSR)
}





















