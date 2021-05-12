#' constructor of EU portfolio object. Type=traditional
#'
#' @inheritParams EUShrinkPortfolio
#' @return an object of class ExUtil_portfolio
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the assets |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean vector estimate for the assets |
#' | W_EU_hat | portfolio weights_sample estimate |
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
#' test <- new_ExUtil_portfolio_traditional(x=x, gamma=gamma)
#' str(test)
#' @export
new_ExUtil_portfolio_traditional <- function(x, gamma){

  if (is.data.frame(x)) x <- as.matrix(x)
  cl <- match.call()
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Portfolio weights
  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfol variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- means %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
  class = c("ExUtil_portfolio"))
}
