
#### Example of input ####
n<-375 # number of realizations
p<-300 # number of assets
b<-rep(1/p,p)
gamma<-5
alph <- 0.05

#set.seed(2)
#x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#### The function ####

set.seed(1)
X <- mvrnorm(n,mu, cov.mat.H0)#<--- simulated matrix of asset returns

# Example of usage
portf_old<- new_ExUtil_portfolio_weights_BDOPS20(x=x, gamma=gamma, b=b)
portf_new <- new_ExUtil_portfolio_weights_BDOPS20_new(x=X, gamma, b, alph)

sum(head(portf_old$means) == head(portf_new$means))
sum(head(portf_old$weights) == head(portf_new$weights))
sum(portf_old$cov_mtrx[1,1:6] == portf_new$cov_mtrx[1,1:6])
sum(head(portf_old$W_EU_hat) == head(portf_new$W_EU_hat))
head(portf_old$alpha) == head(portf_new$alpha)
