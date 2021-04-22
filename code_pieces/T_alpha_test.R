n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

T_alpha <- T_alpha(gamma=gamma, x=x, w_0=b, beta=0.05)
T_alpha


p_value <- 2*(1-pnorm(abs(T_alpha), mean = 0, sd = 1))
p_value


beta <- 0.05

z <- qnorm(p=1-beta/2 , mean = 0, sd = 1)



