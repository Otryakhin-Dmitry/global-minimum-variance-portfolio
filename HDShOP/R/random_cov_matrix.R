# version for tests with arbitrary covariance matrices
SRandCovMtrx <- function(n=3e2, p=2e2, q=20.55){

  cc <- p/n
  ones<-matrix(1,p,1)
  tones<-t(ones)

  #####---Covariance matrix from Wishart distr. with given eigenvalues---#####

  lambda <- rep(NA, p)

  for(i in 1:p){
    lambda[i] <- 0.1*exp(q*cc*(i-1)/p)
  }

  Eigen.matrix <- diag(lambda)

  sq.Eigen.matrix <- sqrt(Eigen.matrix)
  ###---Wishart distribution----####
  Z0 <- matrix(rnorm(p*n*n),p,n*n)
  Wish.matr <- Z0%*%t(Z0)
  ####----spectral decomposition----####
  U <- eigen(Wish.matr)$vectors
  #####---Covariance matrix from Wishart distr. with given eigenvalues---#####
  cov.mat.H0 <- U%*% Eigen.matrix %*%t(U)
  cov.mat.H0
}



#RandCovMtrx <- function(n=3e2, p=2e2, q=20.55){

#  cc <- p/n
  #####---Covariance matrix from Wishart distr. with given eigenvalues---#####

#  eig = 0.1*exp(q*cc*seq(0,1,length=p)) # eig.v[i.eig] was replaced

#  Z  <- array(rnorm(p*p,0,1), dim=c(p,p))
#  QR <- qr(Z)
#  Q  <- qr.Q(QR)
#  R  <- qr.R(QR)
#  D  <- diag(as.vector(diag(R)))
#  D  <-D%*%solve(chol(D%*%D))
#  H  <-Q%*%D
 # E  <-diag(eig)
#  Sigma<-H%*%E%*%t(H)
 # Sigma
#}



# Yarema's new version
#
#' Covariance matrix generator
#'
#' Description: generates a covariance matrix from Wishart distribution with given
#' eigenvalues or with exponentially decreasing eigenvalues. Useful for examples
#' and tests when an arbitrary covariance is needed.
#'
#' This function generates a symmetric positive definite covariance matrix with
#' given eigenvalues. The eigenvalues can be specified either explicitly and by
#' default are generated with exponential decay.
#'
#' @param p dimension of the covariance matrix
#' @param eigenvalues the vector of positive eigenvalues
#' @examples
#'
#' p<-3e2
#' # A non-diagonal covariance matrix
#' Mtrx <- RandCovMtrx(p=p)
#' Mtrx[1:6,1:6]
#' @export
RandCovMtrx <- function(p=2e2, eigenvalues = 0.1*exp(5*seq(0,1,length=p))){

  #####---Covariance matrix from Wishart distr. with given eigenvalues---#####

  if(any(!(eigenvalues > 0))) stop('eigenvalues should be numeric and > 0')
  Eigen.matrix <- diag(eigenvalues)

  ###---Wishart distribution----####
  Z0 <- matrix(rnorm(p*p*p),p,p*p)
  Wish.matr <- Z0%*%t(Z0)
  ####----spectral decomposition----####
  U <- eigen(Wish.matr)$vectors
  #####---Covariance matrix from Wishart distr. with given eigenvalues---#####
  cov.mat.H0 <- U%*% Eigen.matrix %*%t(U)
  cov.mat.H0
}
