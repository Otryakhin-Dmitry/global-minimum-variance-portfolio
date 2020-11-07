

#### constructor of EU portfolio object. The gang, IEEE 2020
new_ExUtil_portfolio <- function(x, gamma, b){

  p <- nrow(x)
  n <- ncol(x)
  cov_mtrx <- Sigma_sample_estimator(x)
  means <- .rowMeans(x, m=p, n=n)

  invSS <- solve(cov_mtrx)
  Ip <- rep.int(1, nrow(x))
  Q_n_hat <- Q_hat_n(x) # this could be optimized

  W_EU_hat <- (invSS %*% Ip)/as.numeric(t(Ip) %*% invSS %*% Ip) +
    Q_n_hat %*% means/gamma

  al <- alpha_hat_star_c(gamma, x=x, b=b)
  weights <- al*W_EU_hat + (1-al)*b
  # methods for covar and means could be mixed in W_EU_hat only, not here nor in BFGSE
  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 weights=weights,
                 alpha=al),
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



#### constructor of EU portfolio object. Type=mean, Mean.type=Bayes-Stein.
new_ExUtil_portfolio_mean_BayesStein <- function(x, gamma, b, mu_0=1){

  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # James-Stein mus and alphas
  alp_JS_hat <- as.numeric((p+2) / (p+2 + n*t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)))
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect

  W_EU_hat <- (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
              Q_n_hat %*% mu_hat_JS / gamma # !!!! gamma instead of alpha^-1 like in Yarema's notes
  al <- alpha_hat_star_c(gamma, x, b)
  W_hat_BFGSE <- al*W_EU_hat + (1-al)*b

  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 mu_hat_JS=mu_hat_JS,
                 alp_JS_hat=alp_JS_hat,
                 W_EU_hat=W_EU_hat,
                 weights=W_hat_BFGSE,
                 alpha=al),
            class = c("ExUtil_portfolio","ExUtil_portfolio_mean_Bayes-Stein")) # add alpha, stand dev, p-value when type=weights
}

# Example
n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
rm(n,p)

test <- new_ExUtil_portfolio_mean_BayesStein(x=x, gamma=gamma, b=b)
str(test)


#### constructor of EU portfolio object. Type=mean, Mean.type=James-Stein.
new_ExUtil_portfolio_mean_JamesStein <- function(x, gamma, b, mu_0=1){

  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # James-Stein mus and alphas
  val <- as.numeric( (p-2) / n / (t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)) )
  alp_JS_hat <- min(1,val)
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect

  W_EU_hat <- (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
    Q_n_hat %*% mu_hat_JS / gamma # !!!! gamma instead of alpha^-1 like in Yarema's notes
  al <- alpha_hat_star_c(gamma, x, b)
  W_hat_BFGSE <- al*W_EU_hat + (1-al)*b

  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 mu_hat_JS=mu_hat_JS,
                 alp_JS_hat=alp_JS_hat,
                 W_EU_hat=W_EU_hat,
                 weights=W_hat_BFGSE,
                 alpha=al),
            class = c("ExUtil_portfolio","ExUtil_portfolio_mean_Bayes-Stein")) # add alpha, stand dev, p-value when type=weights
}

# Example
n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
rm(n,p)

test <- new_ExUtil_portfolio_mean_JamesStein(x=x, gamma=gamma, b=b)
str(test)

