
#' Bayes-Stein shrinkage estimator of the mean vector
#'
#' Bayes-Stein shrinkage estimator of the mean vector as suggested in \insertCite{Jorion1986;textual}{HDShOP}.
#' The estimator is given by \deqn{\hat \mu_{BS} = (1-\beta) \bar x + \beta Y_0 1 \quad ,}
#' where \eqn{\bar x} is the ordinary sample mean, \eqn{\beta} and \eqn{Y_0} are
#' derived using Bayesian approach (see Eq.14 and Eq.17 in \insertCite{Jorion1986;textual}{HDShOP}).
#'
#' @param x a numeric data matrix. Rows represent different variables, columns- observations.
#' @return a numeric vector containing the Bayes-Stein shrinkage estimation of the mean vector
#' @references \insertAllCited{}
#' @examples
#' n <- 7e2 # number of realizations
#' p <- .5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_bs(x=x)
#' @export
mean_bs <- function(x)
{
  p <- nrow(x)
  n <- ncol(x)

  cov_mtrx <- Sigma_sample_estimator(x) * (n-1)/(n-p-2)
  invSS <- solve(cov_mtrx)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  mu_0 <- as.numeric((t(I_vect) %*% invSS %*% means)/(t(I_vect) %*% invSS %*% I_vect))

  # Bayes-Stein mus and alphas
  alp_JS_hat <- as.numeric((p+2) / (p+2 + n*t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)))
  mu_hat_BS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect
  list(means=mu_hat_BS, alpha=alp_JS_hat)
}


#' James-Stein shrinkage estimator of the mean vector
#'
#' James-Stein shrinkage estimator of the mean vector as suggested in \insertCite{Jorion1986;textual}{HDShOP}.
#' The estimator is given by \deqn{\hat \mu_{JS} = (1-\beta) \bar x + \beta Y_0 1 \quad ,}
#' where \eqn{\bar x} is the sample mean vector, \eqn{\beta} is the shrinkage
#' coefficient which minimizes a quadratic loss given by Eq.(11) in \insertCite{Jorion1986;textual}{HDShOP}.
#' \eqn{Y_0} is a prespecified value.
#'
#' @inheritParams mean_bs
#' @param Y_0 a numeric variable. Shrinkage target coefficient.
#' @return a numeric vector containing the James-Stein shrinkage estimator of the mean vector.
#' @references \insertAllCited{}
#' @examples
#' n<-7e2 # number of realizations
#' p<-.5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_js(x=x, Y_0 = 1)
#' @export
mean_js <- function(x, Y_0 = 1)
{
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)

  # James-Stein mus and alphas
  val <- as.numeric( (p-2) / n / (t(means-Y_0*I_vect)%*%invSS%*%(means-Y_0*I_vect)) )
  alp_JS_hat <- min(1,val)
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * Y_0 * I_vect
  list(means=mu_hat_JS, alpha=alp_JS_hat)
}

#' BOP shrinkage estimator
#'
#' Shrinkage estimator of the high-dimensional mean vector as suggested in \insertCite{BOP2019;textual}{HDShOP}.
#' It uses the formula
#' \deqn{\hat \mu_{BOP} = \hat \alpha \bar x + \hat \beta \mu_0 \quad ,} where
#' \eqn{\hat \alpha} and \eqn{\hat \beta} are shrinkage coefficients given by Eq.(6) and Eg.(7)
#' of \insertCite{BOP2019;textual}{HDShOP} that minimize weighted quadratic loss for a given
#' target vector \eqn{\mu_0} (shrinkage target). \eqn{\bar x} stands for the
#' sample mean vector.
#'
#' @param x a p by n matrix or a data frame. Rows represent different variables, columns- observations.
#' @param mu_0 a numeric vector. The target vector used in the construction of the shrinkage estimator.
#' @return a numeric vector containing the shrinkage estimation of the mean vector
#' @references \insertAllCited{}
#' @examples
#' n<-7e2 # number of realizations
#' p<-.5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_bop19(x=x)
#' @export
mean_bop19 <- function(x, mu_0 = rep(1,p))
{
  n <- ncol(x)
  p <- nrow(x)
  samp_mean <- .rowMeans(x=x, m=p, n=n)
  Sigma_n_inv <- solve(Sigma_sample_estimator(x=x))

  alpha_n <- alpha_star_hat_BOP19(n=n, p=p, y_n_aver=samp_mean,
                                  Sigma_n_inv=Sigma_n_inv, mu_0=mu_0)

  beta_n <- beta_star_hat_BOP19(n=n, p=p, alpha_star_hat=alpha_n,
                                y_n_aver_t = t(samp_mean),
                                Sigma_n_inv=Sigma_n_inv, mu_0=mu_0)
  output <-alpha_n *samp_mean + beta_n * mu_0
  list(means=output, alpha_n=alpha_n, beta_n=beta_n)
}
