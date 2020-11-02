
#### Tests for eigen-values

library('MASS')
n<-5e2 # number of realizations
p<-0.4*n # number of assets
w_0 <- rep(1/p,p)
mu <- seq(0.2,-0.2, length.out=p)
q=20.55

Sigma <- SRandCovMtrx(n=n, p=p, mu=mu)
Sigma[1:10,1:10]


cc <- p/n
ones<-matrix(1,p,1)
tones<-t(ones)

lambda <- rep(NA, p)

for(i in 1:p){
lambda[i] <- 0.1*exp(q*cc*(i-1)/p)
}

Eigen.matrix <- diag(lambda)
sq.Eigen.matrix <- sqrt(Eigen.matrix)

lambda
eigen(Sigma)$values[1:10]










