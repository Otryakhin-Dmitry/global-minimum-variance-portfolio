

#### Here we copy Ostap's functions ###########

our.oracle = function(e.mean, t.mean, mu0, icov.matr, params = FALSE){
  t1 = Sys.time()
  # i.cov = ginv(cov.matr, tol = tol.ginv)
  i.cov = icov.matr
  eie = t(e.mean) %*% i.cov %*% e.mean
  ei0 = t(e.mean) %*% i.cov %*% mu0
  ni0 = t(mu0)    %*% i.cov %*% mu0
  eit = t(e.mean) %*% i.cov %*% t.mean
  ti0 = t(t.mean) %*% i.cov %*% mu0
  denominator = eie * ni0 - ei0^2
  alpha_n = as.double((eit %*% ni0 - ti0 %*% ei0) / denominator)
  beta_n  = as.double((eie %*% ti0 - ei0 %*% eit) / denominator)
  if(!params){
    c(as.vector(alpha_n * e.mean + beta_n * mu0), difftime(Sys.time(), t1, units = "secs"))
  }else{
    c(alpha = alpha_n, beta = beta_n)
  }
}

####

our.bonafide = function(X, e.mean, ie.cov, mu0, params = FALSE){ ## gooood!
  t1 = Sys.time()
  p     = dim(X)[2]
  n     = dim(X)[1]
  cc    = p / n
  i.cov = ie.cov
  # i.cov = ginv(e.cov, tol = tol.ginv)

  q00 = t(mu0) %*% i.cov %*% mu0
  qnn = t(e.mean) %*% i.cov %*% e.mean
  q0n = t(e.mean) %*% i.cov %*% mu0

  # gcc  = if(cc < 1){cc / (1 - cc)}else{1 / (cc - 1)}
  gcc  = if(cc < 1){cc / (1 - cc)}else{1 / (cc - 1)}
  alpha_hat = as.double(((qnn - gcc) * q00 - q0n^2) / (qnn * q00 - q0n^2))
  beta_hat = as.double((1 - alpha_hat) * (q0n / q00))
  if(!params){
    c(as.vector(alpha_hat * e.mean + beta_hat * mu0), difftime(Sys.time(), t1, units = "secs"))
  }else{
    c(alpha = alpha_hat, beta = beta_hat)
  }
}

#### Parameter setting ####

c <- 0.3
n <- 5e2
p <- c*n

#### Computing the results ####

mu_0<-rep(5, p)
mu_n<-rep(10, p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 5
invSigma <- solve(Sigma)


#### Simulation and computation is here ####
dv <- MASS::mvrnorm(n = n, mu = mu_n, Sigma = Sigma)
X <- matrix(data=dv, nrow=p, ncol=n, byrow = FALSE)

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(cov(X))

our.oracle(e.mean=y_n_aver, t.mean=mu_n, mu0=mu_0, icov.matr=invSigma, params = TRUE)

our.bonafide(X=X, e.mean=y_n_aver, ie.cov=Sigma_n_inv, mu0=mu_0, params = TRUE)


####
dv <- MASS::mvrnorm(n = n, mu = mu_n, Sigma = Sigma)
X <- matrix(data=dv, nrow=p, ncol=n, byrow = FALSE)

y_n_aver <- rowMeans(X)


alpha_star_n_BOP19(y_n_aver=y_n_aver, Sigma_n_inv=invSigma, mu_n=mu_n, mu_0=mu_0)
al<-alpha_star_BOP19(c=c, mu_n=mu_n, Sigma_n_inv=invSigma, mu_0=mu_0)
al
# al<-alpha_star_hat_BOP19(n=n, p=p, y_n_aver=y_n_aver, Sigma_n_inv=Sigma, mu_0=mu_0)
# al


beta_star_n_BOP19(y_n_aver=y_n_aver, Sigma_n_inv=invSigma, mu_n=mu_n, mu_0=mu_0)
beta_star_BOP19(alpha_star=al, mu_n_t=t(mu_n), Sigma_n_inv=invSigma, mu_0=mu_0)
# beta_star_hat_BOP19(n=n, p=p, alpha_star_hat=al,
                    # y_n_aver_t=t(y_n_aver), Sigma_n_inv=Sigma, mu_0=mu_0)


our.oracle(e.mean=y_n_aver, t.mean=mu_n, mu0=mu_0, icov.matr=Sigma, params = TRUE)
our.bonafide(X=X, e.mean=y_n_aver, ie.cov=Sigma, mu0=mu_0, params = TRUE)

