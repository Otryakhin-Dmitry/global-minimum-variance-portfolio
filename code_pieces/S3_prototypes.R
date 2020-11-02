

# constructor of EU portfolio object
new_ExUtil_portfolio <- function(x, gamma, b){

  cov_mtrx <- Sigma_sample_estimator(x)
  means <- .rowMeans(x, m=p, n=n) # could be computed via James-Stein
  weights <- W_hat_BFGSE(x=x, gamma=gamma, b=b)
 # W_EU_hat to return as well
 # metthods for covar and means could be mixed in W_EU_hat only, not here nor in BFGSE
  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 weights=weights),
            class = "ExUtil_portfolio") # add alpha, stand dev, p-value when type=weights
}



# Example

n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

test <- new_ExUtil_portfolio(x=x, gamma=gamma, b=b)
str(test)




