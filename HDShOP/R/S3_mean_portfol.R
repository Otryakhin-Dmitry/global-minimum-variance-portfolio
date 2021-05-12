
#' constructor of EU portfolio object. Type=mean, Mean.type=Bayes-Stein.
#'
#' @inheritParams EUShrinkPortfolio
#' @param mu_0 a numeric value. The scaling of the target for shrinkage of the mean vector.
#' @references \insertRef{Jorion1986}{HDShOP}
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the assets |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | mean vector estimate for the assets obtained through the Bayes-Stein method |
#' | alp_BS_hat | shrinkage intensity for the mean vector |
#' | W_EU_hat | portfolio weights estimate computed via the above means and cov_mtrx |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_ExUtil_portfolio_mean_BayesStein(x=x, gamma=gamma, mu_0 = 0)
#' str(test)
#' @export
new_ExUtil_portfolio_mean_BayesStein <- function(x, gamma, mu_0=0){

  cl <- match.call()
  if (is.data.frame(x)) x <- as.matrix(x)
  p <- nrow(x)
  n <- ncol(x)

  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Bayes-Stein mus and alphas
  alp_BS_hat <- as.numeric((p+2) / (p+2 + n*t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)))
  mu_hat_BS <- (1-alp_BS_hat) * means + alp_BS_hat * mu_0 * I_vect

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
    Q_n_hat %*% mu_hat_BS / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- mu_hat_BS %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=mu_hat_BS,
                 alp_BS_hat=alp_BS_hat,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("ExUtil_portfolio_mean_Bayes-Stein", "ExUtil_portfolio")) # add alpha, stand dev, p-value when type=weights
}


#' constructor of EU portfolio object. Type=mean, Mean.type=James-Stein.
#'
#' @inheritParams new_ExUtil_portfolio_mean_BayesStein
#' @references \insertRef{Jorion1986}{HDShOP}
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the assets |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | mean vector estimate for the assets obtained through the James-Stein method |
#' | alp_JS_hat | shrinkage intensity for the mean vector |
#' | W_EU_hat | portfolio weights estimate computed via the above means and cov_mtrx |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_ExUtil_portfolio_mean_JamesStein(x=x, gamma=gamma, mu_0=0)
#' str(test)
#' @export
new_ExUtil_portfolio_mean_JamesStein <- function(x, gamma, mu_0=0){

  cl <- match.call()
  if (is.data.frame(x)) x <- as.matrix(x)

  p <- nrow(x)
  n <- ncol(x)

  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # James-Stein mus and alphas
  val <- as.numeric( (p-2) / n / (t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)) )
  alp_JS_hat <- min(1,val)
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
    Q_n_hat %*% mu_hat_JS / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- mu_hat_JS %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=mu_hat_JS,
                 alp_JS_hat=alp_JS_hat,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("ExUtil_portfolio_mean_James-Stein", "ExUtil_portfolio")) # add alpha, stand dev, p-value when type=weights
}


#' constructor of EU portfolio object. Type=mean, Mean.type=BOP.
#'
#' @inheritParams EUShrinkPortfolio
#' @param mu_0 a numeric vector. The target for shrinkage of the mean vector of asset returns.
#' @references \insertRef{BOP2019}{HDShOP}
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the assets |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | mean vector estimate for the assets obtained through the method by Bodnar et al |
#' | alpha_star_hat | shrinkage intensity for the mean vector |
#' | beta_star_hat | shrinkage intensity for the target vector |
#' | W_EU_hat | portfolio weights estimate computed via the above means and cov_mtrx |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#' mu_0 <- rep(1,p)
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_ExUtil_portfolio_mean_BOP19(x=x, gamma=gamma, mu_0=mu_0)
#' str(test)
#' @export
new_ExUtil_portfolio_mean_BOP19 <- function(x, gamma, mu_0=mu_0){

  cl <- match.call()
  if (is.data.frame(x)) x <- as.matrix(x)

  p <- nrow(x)
  n <- ncol(x)

  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Bodnar-Okhrin-Parolya mus, alpha and beta
  alpha_star_hat <- alpha_star_hat_BOP19(n=n, p=p, y_n_aver=means, Sigma_n_inv=invSS, mu_0=mu_0)
  beta_star_hat <- beta_star_hat_BOP19(n=n, p=p, alpha_star_hat=alpha_star_hat, y_n_aver_t=t(means), Sigma_n_inv=invSS, mu_0=mu_0)
  mu_hat_BOP <- alpha_star_hat * means + beta_star_hat * mu_0

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% mu_hat_BOP / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- mu_hat_BOP %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=mu_hat_BOP,
                 alpha_star_hat=alpha_star_hat,
                 beta_star_hat=beta_star_hat,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
  class = c("ExUtil_portfolio_mean_BOP19", "ExUtil_portfolio")) # add alpha, stand dev, p-value when type=weights
}
