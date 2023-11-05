
# Yarema's new version
#
#' Covariance matrix generator
#'
#' Generates a covariance matrix from Wishart distribution with given
#' eigenvalues or with exponentially decreasing eigenvalues. Useful for examples
#' and tests when an arbitrary covariance matrix is needed.
#'
#' This function generates a symmetric positive definite covariance matrix with
#' given eigenvalues. The eigenvalues can be specified explicitly. Or, by
#' default, they are generated with exponential decay.
#'
#' @param p dimension of the covariance matrix
#' @param eigenvalues the vector of positive eigenvalues
#' @return covariance matrix
#' @examples
#'
#' p<-1e1
#' # A non-diagonal covariance matrix
#' Mtrx <- RandCovMtrx(p=p)
#' Mtrx
#' @export
RandCovMtrx <- function(p=2e2, eigenvalues = 0.1*exp(5*seq_len(p)/p)){

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
