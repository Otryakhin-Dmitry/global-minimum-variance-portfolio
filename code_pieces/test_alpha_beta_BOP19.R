c <- 0.5
n <- 2e3
p <- c*n

mu_0<-rep(1/sqrt(p), p)
mu_n<-rep(1/sqrt(p), p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- solve(Sigma)

X <- matrix(data=rnorm(n*p), nrow=p, ncol=n) # replace with Y matrix. X-

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)

alpha_star_n_BOP19(y_n_aver, Sigma_n_inv, mu_n, mu_0)
alpha_star_hat_BOP19(n, p, y_n_aver, Sigma_n_inv, mu_0)

s <- s_BOP19(mu_0, Sigma_n_inv, mu_n)
sigma_s_square <- sigma_s_square_BOP19(c=c, s=s)
R_BOP19 <- R_BOP19(mu_0, Sigma_n_inv, mu_n)
Omega_BOP19(c, sigma_s_square, s, R_BOP19)


#### Diagnostics  #################

library(EstimDiagnostics)
c<-0.5
Nmc=400
s<-c(1e2,1e3)

Inference<-function(s){

  n <- s
  p <- c*n

  mu_0<-rep(1, p)
  mu_n<-rep(0, p)
  Sigma=matrix(0, p, p)
  diag(Sigma) <- 1
  invSigma <- solve(Sigma)

  X <- matrix(data=rnorm(p*n), nrow=p, ncol=n)

  y_n_aver <- rowMeans(X)
  Sigma_n_inv <- Sigma_sample_estimator(X)

  alpha_star_n <- alpha_star_n_BOP19(y_n_aver=y_n_aver, Sigma_n_inv=invSigma,
                                     mu_n = mu_n, mu_0 = mu_0)
  beta_star_n <-beta_star_n_BOP19(y_n_aver=y_n_aver, Sigma_n_inv=invSigma,
                                        mu_n = mu_n, mu_0 = mu_0)
  alpha_star_hat <- alpha_star_hat_BOP19(n=n, p=p, y_n_aver=y_n_aver,
                                         Sigma_n_inv=Sigma_n_inv, mu_0=mu_0)
  beta_star_hat <- beta_star_hat_BOP19(n=n, p=p, alpha_star_hat=alpha_star_hat,
                                       y_n_aver_t=t(y_n_aver), Sigma_n_inv=Sigma_n_inv,
                                       mu_0=mu_0)
  alph_dif <- alpha_star_n - alpha_star_hat
  beta_dif <- beta_star_n - beta_star_hat

  list(alph_dif=alph_dif, beta_dif=beta_dif)
}

data <- Estim_diagnost(Nmc, s, Inference)
estims_qqplot(data)
estims_boxplot(data)









