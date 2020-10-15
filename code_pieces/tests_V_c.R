
###############################

n<-8e2 # number of realizations
p<-0.8*n # number of assets
w_0 <- rep(1/p,p)
gamma<-5e4

###############################


x <-matrix(data = rnorm(n=n*p), nrow = p, ncol = n)

V_hat_c(x)#*(1-p/n)
V_hat_GMV(x)











