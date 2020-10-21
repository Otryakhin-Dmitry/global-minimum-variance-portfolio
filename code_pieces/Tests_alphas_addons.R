

library('MASS')

n<-5e2 # number of realizations
p<-0.4*n # number of assets
w_0 <- rep(1/p,p)
mu <- seq(0.2,-0.2, length.out=p)

Sigma <- SRandCovMtrx(n=n, p=p, mu=mu)


vect_as <- sqrt(n)*
  replicate(n=3e2, {
    x <- t(mvrnorm(n=n, mu=mu, Sigma=Sigma))
    alpha_hat_star_c_GMV(x, b=w_0) -
      alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)
  })


(sd(vect_as))^2
SD_alpha_simple(Sigma=Sigma, b=w_0, mu=mu, n=n)




