#' Traditional mean-variance portfolio
#'
#' Mean-variance portfolios with the traditional (sample) estimators for the mean
#' vector and the covariance matrix of asset returns. For more details of the method,
#' see \code{\link{MVShrinkPortfolio}}. new_MV_portfolio_traditional is for the
#' case p<n, while new_MV_portfolio_traditional_pgn is for p>n, where p is the
#' number of assets and n is the number of observations.
#'
#' @inheritParams MVShrinkPortfolio
#' @return an object of class MeanVar_portfolio
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of asset returns |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean estimator of the asset returns |
#' | W_mv_hat | sample estimator of portfolio weights |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | expected portfolio return |
#' | Sharpe | portfolio Sharpe ratio |
#' @md
#'
#' @examples
#' n <- 3e2 # number of realizations
#' p <- .5*n # number of assets
#' gamma <- 1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_MV_portfolio_traditional(x=x, gamma=gamma)
#' str(test)
#' @export
new_MV_portfolio_traditional <- function(x, gamma){

  if (is.data.frame(x)) x <- as.matrix(x)
  cl <- match.call()
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/
                     as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Portfolio weights
  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfolio variance
  Port_Var <- as.numeric(t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat)
  Port_mean_return <- as.numeric(means %*% W_EU_hat)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 weights=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
  class = c("MeanVar_portfolio"))
}


#' @rdname new_MV_portfolio_traditional
#' @export
new_MV_portfolio_traditional_pgn <- function(x, gamma){

  if (is.data.frame(x)) x <- as.matrix(x)
  cl <- match.call()

  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- MASS::ginv(cov_mtrx)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/
                     as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Portfolio weights
  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  # Sharpe, mean return and portfolio variance
  Port_Var <- as.numeric(t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat)
  Port_mean_return <- as.numeric(means %*% W_EU_hat)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 weights=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("MeanVar_portfolio"))
}
