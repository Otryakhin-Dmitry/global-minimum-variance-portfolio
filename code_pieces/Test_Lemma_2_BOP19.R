

#### Examples and tests ####
n<-3e2
c<-0.3
p<-c*n

mu_0<-rep(1/sqrt(p), p)
mu_n<-rep(1/sqrt(p), p)

# Random matrix
#Sigma <- RandCovMtrx(n=n, p=p, q=20.55, mu=mu_n)

# I matrix
Sigma <- matrix(0, p, p)
diag(Sigma) <- 1

sqrM <- expm::sqrtm(Sigma)

X <- t(MASS::mvrnorm(n=n, mu=mu_n, Sigma=Sigma))
invS <- solve(Sigma_sample_estimator(X))
x_bar_n <- rowMeans(X) - mu_n

t(x_bar_n) %*% sqrM %*% invS %*% sqrM %*%  x_bar_n
c/(1-c)

1/p * t(x_bar_n) %*% Sigma %*% invS %*% mu_n

1/p * abs(t(mu_n) %*% invS %*% mu_n - 1/(1-c) * t(mu_n) %*% Sigma %*% mu_n)


