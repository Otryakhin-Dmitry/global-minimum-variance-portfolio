c <- 0.3
n <- 2e3
p <- c*n

mu_0<-rep(10, p)
mu_n<-rep(0, p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- Sigma #solve(Sigma)


X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)

alpha_star_n_BOP19(y_n_aver, Sigma_n_inv, mu_n, mu_0)
alpha_star_hat_BOP19(n, p, y_n_aver, Sigma_n_inv, mu_0)
