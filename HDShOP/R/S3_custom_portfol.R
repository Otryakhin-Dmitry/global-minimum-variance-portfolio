
#' A constructor for class MeanVar_portfolio
#'
#' A light-weight constructor of objects of S3 class MeanVar_portfolio. This  function is for development purposes.
#' A helper function equipped with error messages and allowing more flexible input is \code{\link{MeanVar_portfolio}}.
#'
#' @param mean_vec mean vector of asset returns
#' @param cov_mtrx the covariance matrix of asset returns
#' @inheritParams MVShrinkPortfolio
#' @return Mean-variance portfolio in the form of object of S3 class MeanVar_portfolio.
#'
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' # Simple MV portfolio
#' cov_mtrx <- Sigma_sample_estimator(x)
#' means <- rowMeans(x)
#'
#' cust_port_simp <- new_MeanVar_portfolio(mean_vec=means, cov_mtrx=cov_mtrx, gamma=2)
#' str(cust_port_simp)
#'
#' # Portfolio with Bayes-Stein shrunk means
#' # and a Ledoit and Wolf estimator for covariance matrix
#' TM <- matrix(0, p, p)
#' diag(TM) <- 1
#' cov_mtrx <- CovarEstim(x, type="LW20", TM=TM)
#' means <- rowMeans(x)
#'
#' cust_port_BS_LW <- new_MeanVar_portfolio(mean_vec=means, cov_mtrx=cov_mtrx, gamma=2)
#' str(cust_port_BS_LW)
#' @export
new_MeanVar_portfolio <- function(mean_vec, cov_mtrx, gamma){

  cl <- match.call()

  inv_cov_mtrx <- solve(cov_mtrx)
  p <- nrow(inv_cov_mtrx)
  I_vect <- rep(1, times=p)

  Q_n_hat <- inv_cov_mtrx - (inv_cov_mtrx %*% I_vect %*% t(I_vect) %*% inv_cov_mtrx)/as.numeric(t(I_vect) %*% inv_cov_mtrx %*% I_vect)

  W_EU_hat <- as.vector(
    (inv_cov_mtrx %*% I_vect)/as.numeric(t(I_vect) %*% inv_cov_mtrx %*% I_vect) +
      Q_n_hat %*% mean_vec / gamma,
    mode = 'numeric')

  Port_Var <- as.numeric(t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat)
  Port_mean_return <- as.numeric(mean_vec %*% W_EU_hat)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=inv_cov_mtrx,
                 means=mean_vec,
                 weights=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = "MeanVar_portfolio")
}



#' A validator for objects of class MeanVar_portfolio
#' @param w Object of class MeanVar_portfolio.
#' @return If the object passes all the checks, then w itself is returned,
#' otherwise an error is thrown.
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' # Simple MV portfolio
#' cov_mtrx <- Sigma_sample_estimator(x)
#' means <- rowMeans(x)
#'
#' cust_port_simp <- new_MeanVar_portfolio(mean_vec=means, cov_mtrx=cov_mtrx, gamma=2)
#' str(validate_MeanVar_portfolio(cust_port_simp))
#' @export
validate_MeanVar_portfolio <- function(w) {

  values <- unclass(w)

  if (is.null(values$cov_mtrx))  stop("a covariance matrix is missing", call. = FALSE)
  if (is.null(values$inv_cov_mtrx))  stop("an inverse covariance matrix is missing", call. = FALSE)
  if (is.null(values$means))  stop("a mean vector is missing", call. = FALSE)
  if (is.null(values$weights))  stop("a vector of weights is missing", call. = FALSE)

  if (is.null(values$Port_Var))  stop("a portfolio variance is missing", call. = FALSE)
  if (is.null(values$Port_mean_return))  stop("a portfolio mean return is missing", call. = FALSE)
  if (is.null(values$Sharpe))  stop("a Sharpe ratio is missing", call. = FALSE)

  if (!is.vector(values$means))  stop("means is not a vector", call. = FALSE)
  if (!is.vector(values$weights))  stop("weights is not a vector", call. = FALSE)
  if (!identical(class(values$inv_cov_mtrx), c("matrix", "array"))){
    stop("inv_cov_mtrx is not a matrix", call. = FALSE)
  }

  if (length(values$means)!=length(values$weights) | nrow(values$inv_cov_mtrx)!=length(values$weights)) {
    stop(
      "lenghts of the mean vector and the weights must equal the row number of the covariance matrix",
      call. = FALSE
    )
  }

  w
}

#' A helper function for MeanVar_portfolio
#'
#' A user-friendly function making mean-variance portfolios for assets with customly
#' computed covariance matrix and mean returns.
#' The weights are computed in accordance with the formula
#' \deqn{\hat w_{MV} = \frac{\hat{\Sigma}^{-1} 1}{1' \hat{\Sigma}^{-1} 1} + \gamma^{-1} \hat Q \hat{\mu} \quad ,}
#' where \eqn{\hat{\Sigma}} is an estimator for the covariance matrix, \eqn{\hat{\mu}} is an estimator for the mean vector,
#' \eqn{\gamma} is the coefficient of risk aversion, and
#' \eqn{\hat Q} is given by
#' \deqn{\hat Q = \hat{\Sigma}^{-1} - \frac{\hat{\Sigma}^{-1} 1 1' \hat{\Sigma}^{-1}}{1' \hat{\Sigma}^{-1} 1} .}
#' The computation is made by \code{\link{new_MeanVar_portfolio}} and the result
#' is validated by 
#' \code{\link{validate_MeanVar_portfolio}}.
#'
#' @param mean_vec mean vector of asset returns provided in the form of a vector or a list.
#' @param cov_mtrx the covariance matrix of asset returns. It could be a matrix or a data frame.
#' @inheritParams MVShrinkPortfolio
#' @return Mean-variance portfolio in the form of object of S3 class MeanVar_portfolio.
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' # Simple MV portfolio
#' cov_mtrx <- Sigma_sample_estimator(x)
#' means <- rowMeans(x)
#'
#' cust_port_simp <- MeanVar_portfolio(mean_vec=means, cov_mtrx=cov_mtrx, gamma=2)
#' str(cust_port_simp)
#' @export
MeanVar_portfolio <- function(mean_vec, cov_mtrx, gamma){

  cov_mtrx <- as.matrix(cov_mtrx)
  mean_vec<-unlist(mean_vec)

  xx <- new_MeanVar_portfolio(mean_vec=mean_vec,
                              cov_mtrx=cov_mtrx,
                              gamma=gamma)

  validate_MeanVar_portfolio(xx)
}


# Summary method for MeanVar_portfolio
#' @export
summary.MeanVar_portfolio <- function(object, ...){

  list(call=object$call,
       Port_Var=object$Port_Var,
       Port_mean_return=object$Port_mean_return,
       Sharpe=object$Sharpe
       )
}

