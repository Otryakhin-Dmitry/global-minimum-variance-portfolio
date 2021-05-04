#' A constructor for class ExUtil_portfolio
#'
#' A light-weight constructor of objects of S3 class ExUtil_portfolio. This  function is for development purposes.
#' A helper function equipped with error messages and allowing more flexible input is \code{\link{ExUtil_portfolio_custom}}.
#'
#' @param mean_vec mean vector of asset returns
#' @param inv_cov_mtrx the inverse covariance matrix of asset returns
#' @inheritParams EUShrinkPortfolio
#' @return Object of S3 class ExUtil_portfolio.
#'
#' This class is designed to describe Expected Utility portfolios. It comprises portfolio weights
#' W_EU_hat, mean vector means and an inverse covariance inv_cov_mtrx. W_EU_hat and means both
#' have the form of a numeric vector, while inv_cov_mtrx is a matrix. The direct covariance matrix
#' is not included in order to cover cases where only the inverse one exists.
#'
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

  structure(list(inv_cov_mtrx=inv_cov_mtrx,
                 means=mean_vec,
                 W_EU_hat=W_EU_hat
                 ),
            class = "ExUtil_portfolio")
}

#' A validator for ExUtil_portfolio
#' @param x Object of class ExUtil_portfolio.
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
#' str(validate_ExUtil_portfolio(cust_port_simp))
#' @export
validate_ExUtil_portfolio <- function(x) {

  values <- unclass(x)

  if (is.null(values$inv_cov_mtrx))  stop("an inverse covariance matrix is missing", call. = FALSE)
  if (is.null(values$means))  stop("a mean vector is missing", call. = FALSE)
  if (is.null(values$W_EU_hat))  stop("a vector of weights is missing", call. = FALSE)

  if (!is.vector(values$means))  stop("means is not a vector", call. = FALSE)
  if (!identical(class(values$inv_cov_mtrx), c("matrix", "array"))){
    stop("inv_cov_mtrx is not a matrix", call. = FALSE)
  }

  if (length(values$means)!=length(values$W_EU_hat) | nrow(values$inv_cov_mtrx)!=length(values$W_EU_hat)) {
    stop(
      "lenghts of the mean vector and the weights must equal the row number of the covariance matrix",
      call. = FALSE
    )
  }

  x
}

#' A helper function for ExUtil_portfolio
#' @param mean_vec mean vector or list of asset returns.
#' @param inv_cov_mtrx the inverse covariance matrix of asset returns. Could be a matrix or a data frame.
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
#' cust_port_simp <- ExUtil_portfolio_custom(mean_vec=means, inv_cov_mtrx=invSS, gamma=2)
#' str(cust_port_simp)
#' @export
ExUtil_portfolio_custom <- function(mean_vec, inv_cov_mtrx, gamma){

  inv_cov_mtrx <- as.matrix(inv_cov_mtrx)
  mean_vec<-unlist(mean_vec)

  xx <- new_ExUtil_portfolio_custom(mean_vec=mean_vec,
                                    inv_cov_mtrx=inv_cov_mtrx,
                                    gamma=gamma)
  validate_ExUtil_portfolio(xx)
}


# Summary method for ExUtil_portfolio
#' @export
summary.ExUtil_portfolio <- function(object){

  list(call=object$call)
}

#' @export
summary.ExUtil_portfolio_weights_BDOPS20 <- function(object){

  list(call=object$call,
       Port_Var=object$Port_Var,
       Port_mean_return=object$Port_mean_return,
       Sharpe=object$Sharpe
      )
}

#' @export
summary.GMV_portfolio_weights_BDPS19 <- function(object){

  list(call=object$call,
       Port_Var=object$Port_Var,
       Port_mean_return=object$Port_mean_return,
       Sharpe=object$Sharpe
  )
}
