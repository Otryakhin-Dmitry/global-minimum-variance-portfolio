

#### Parameter setting
n<-1e3 # number of realizations
p<-0.8*n # number of assets
w_0 <- rep(1/p,p)
# gamma<-5e4


#### alphas EU work when gamma=Inf
x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous


Sigma <- matrix(0, p, p)
diag(Sigma) <- 1
alpha_star(gamma=Inf, mu=rep(0,p), Sigma=Sigma, c=p/n, b=w_0) # must be computable
alpha_star_GMV(mu=rep(0,p), Sigma=Sigma, c=p/n, b=w_0) # must be equal to the previous


#### Equivalence of alphas for EU and GMV portfolios


#### SD of simple alphas

vect_as <- sqrt(n)*
  replicate(n=1e2, {
    x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
    alpha_hat_star_c_GMV(x, b=w_0) -
      alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)
  })


SD_alpha_simple <- function(Sigma, b, mu, n){

  c <- nrow(Sigma)/n
  V_b <- V_b(Sigma, b)
  V_GMV <- V_GMV(Sigma)
  Lb <- V_b/V_GMV - 1
  R_b <- R_b(mu, b)

  numer <- 2*(1-c)*c^2*(Lb+1)
  denom <- ((1-c)*R_b+c)^4
  multip<- (2-c)*Lb +c

  numer / denom * multip
}

sd(vect_as)
(SD_alpha_simple(Sigma=Sigma, b=w_0, mu=rep(0,p), n=n))^2

abs(
  sd(vect_as*sqrt(n)) -
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






