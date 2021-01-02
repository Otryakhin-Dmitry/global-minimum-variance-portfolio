c <- 0.7
n <- 5e2
p <- c*n

mu_0<-rep(0, p); mu_0[1:4]<-c(100,50,4,1)
mu_n<-rep(0, p); mu_n[1:4]<-c(-10,50,-4,1)
#mu_n<-rep(1/sqrt(p), p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 10
invSigma <- Sigma #solve(Sigma)


# X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)
X <- t(MASS::mvrnorm(n=n, mu=mu_n, Sigma=Sigma))

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)

alpha_star_n_BOP19(y_n_aver, Sigma_n_inv, mu_n, mu_0)
alpha_star_hat_BOP19(n, p, y_n_aver, Sigma_n_inv, mu_0)
alpha_star_BOP19(c=c, mu_n, Sigma_n_inv, mu_0)
