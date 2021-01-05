
#' Custom EU portfolio constructors
#'
#' @param mean_vec mean vector of asset returns
#' @param inv_cov_mtrx the inverse covariance matrix of asset returns
#' @inheritParams EUShrinkPortfolio
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' # Simple EU portfolio
#' cov_mtrx <- Sigma_sample_estimator(x)
#' invSS <- solve(cov_mtrx)
#' means <- rowMeans(x)
#'
#' cust_port_simp <- new_ExUtil_portfolio_custom(mean_vec=means, inv_cov_mtrx=invSS, gamma=2)
#' str(cust_port_simp)
#'
#' # Portfolio with Bayes-Stein shrinked means
#' # and a Ledoit and Wolf estimator for covariance matrix
#' TM <- matrix(0, p, p)
#' diag(TM) <- 1
#' cov_mtrx <- CovarEstim(x, type="LW20", TM=TM)
#' invSS <- solve(cov_mtrx)
#' means <- rowMeans(x)
#'
#' cust_port_BS_LW <- new_ExUtil_portfolio_custom(mean_vec=means, inv_cov_mtrx=invSS, gamma=2)
#' str(cust_port_BS_LW)
#' @export
new_ExUtil_portfolio_custom <- function(mean_vec, inv_cov_mtrx, gamma){

  p <- nrow(inv_cov_mtrx)
  I_vect <- rep(1, times=p)

  Q_n_hat <- inv_cov_mtrx - (inv_cov_mtrx %*% I_vect %*% t(I_vect) %*% inv_cov_mtrx)/as.numeric(t(I_vect) %*% inv_cov_mtrx %*% I_vect)

  W_EU_hat <- as.vector(
    (inv_cov_mtrx %*% I_vect)/as.numeric(t(I_vect) %*% inv_cov_mtrx %*% I_vect) +
       Q_n_hat %*% mean_vec / gamma,
    mode = 'numeric')

  structure(list(W_EU_hat=W_EU_hat),
            class = "ExUtil_portfolio_custom")
}



