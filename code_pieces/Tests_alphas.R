library('MASS')

####

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
  replicate(n=3e2, {
    x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)
    alpha_hat_star_c(gamma=Inf, x=x, b=w_0) -
      alpha_star(gamma=Inf, mu=rep(0,p), Sigma=Sigma, b=w_0, c=p/n)
  })

sd(vect_as)
sd(vect_gen)
SD_alpha_simple(Sigma=Sigma, b=w_0, mu=rep(0,p), n=n)






