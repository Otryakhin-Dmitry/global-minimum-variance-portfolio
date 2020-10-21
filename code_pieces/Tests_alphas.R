library('MASS')

#### Parameter setting
n<-1e3 # number of realizations
p<-0.8*n # number of assets
w_0 <- rep(1/p,p)


#### alphas EU work when gamma=Inf
x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous

# x with dependent components
n<-3e2 # number of realizations
p<-0.3*n # number of assets
w_0 <- rep(1/p,p)
mu <- seq(0.2,-0.2, length.out=p)

cov_mat <- SRandCovMtrx(n=n, p=p, mu=mu)
x <- t(mvrnorm(n,mu, cov_mat))

alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous


Sigma <- matrix(0, p, p)
diag(Sigma) <- 1
alpha_star(gamma=Inf, mu=rep(0,p), Sigma=Sigma, c=p/n, b=w_0) # must be computable
alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0) # must be equal to the previous

#### Equivalence of alphas for EU and GMV portfolios

n<-7e2 # number of realizations
p<-0.4*n # number of assets
w_0 <- rep(1/p,p)

Sigma <- matrix(0, p, p)
diag(Sigma) <- 1

#### SD of simple alphas

vect_as <- sqrt(n)*
  replicate(n=3e2, {
    x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
    alpha_hat_star_c_GMV(x, b=w_0) -
      alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)
  })


sd(vect_as)
SD_alpha_simple(Sigma=Sigma, b=w_0, mu=rep(0,p), n=n)

abs(
  sd(vect_as) -
  SD_alpha_simple(Sigma=Sigma, b=w_0, mu=rep(0,p), n=n)
)


#### Alpha simple computed through the general one ####
vect_gen <- sqrt(n)*
  replicate(n=5e2, {
    x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
    alpha_hat_star_c(gamma=Inf, x=x, b=w_0) -
      alpha_star(gamma=Inf, mu=rep(0,p), Sigma=Sigma, b=w_0, c=p/n)
  })

sd(vect_as)
sd(vect_gen)
SD_alpha_simple(Sigma=Sigma, b=w_0, mu=rep(0,p), n=n)






