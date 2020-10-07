
#### Sigma simple estimator (5) in EU tests paper
#' Simple covariance estimator
#'
#' Some description
#'
#' @param x numeric matrix in which columns are independent realizations of asset returns
#'
#' @examples
#' p<-5 # number of assets
#' n<-1e1 # number of realizations
#'
#' x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' Sigma_sample_estimator(x)
Sigma_sample_estimator <- function(x) {

  a <- rowMeans(x, na.rm = TRUE)
  a_x_size <- matrix(rep(a,ncol(x)),nrow=nrow(x), ncol=ncol(x))
  (x-a_x_size) %*% t(x-a_x_size)/(ncol(x)-1)
}


#### Q estimator (page 3, IEEE)
#' Q_n_hat
#'
#' Q estimator from (6) of EU_IEEE.
#'
#' @inheritParams Sigma_sample_estimator
#'
Q_hat_n <- function(x){

  SS <- Sigma_sample_estimator(x)
  invSS <- solve(SS)
  Ip <- rep.int(1, nrow(x))

  invSS - (invSS %*% Ip %*% t(Ip) %*% invSS)/as.numeric(t(Ip) %*% invSS %*% Ip)
}


#### W_EU estimator (page 3, IEEE)
#' W_EU estimator
#'
#' EU portfolio weights estimator from (6) of EU_IEEE.
#'
#' @inheritParams Sigma_sample_estimator
#' @param gamma positive numeric. Investors attitude towards risk
#'
W_EU_hat <- function(x, gamma){

  SS <- Sigma_sample_estimator(x)
  invSS <- solve(SS)
  Ip <- rep.int(1, nrow(x))
  Q_n_hat <- Q_hat_n(x)

  (invSS %*% Ip)/as.numeric(t(Ip) %*% invSS %*% Ip) +
  Q_n_hat %*% rowMeans(x, na.rm = TRUE)/gamma
}




















