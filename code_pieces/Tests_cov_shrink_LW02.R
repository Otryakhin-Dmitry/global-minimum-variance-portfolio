



#### Examples and tests ####
n<-5e2
c<-0.2
p<-c*n
mu <- rep(0, p)
Sigma <- RandCovMtrx(n=n, p=p, q=20.55, mu=mu)

X <- MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma)
Sigma_shr <- nonlin_shrinkLW(X)

#Sigma[1:10,1:10]
#Sigma_shr[1:10,1:10]
#abs((Sigma_shr[1:10,1:10] - Sigma[1:10,1:10])/Sigma[1:10,1:10])
sum(abs((Sigma_shr[1:5,1:5] - Sigma[1:5,1:5])/Sigma[1:5,1:5]))


