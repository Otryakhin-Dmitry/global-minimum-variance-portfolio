# version for tests with arbitrary covariance matrices
SRandCovMtrx <- function(n=3e2, p=2e2, q=20.55, mu=seq(0.2,-0.2, length.out=p)){

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

