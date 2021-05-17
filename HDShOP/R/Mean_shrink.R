
#' Bayes-Stein shrinkage estimator of the mean vector
#'
#' Bayes-Stein shrinkage estimator of the mean vector as suggested in Jorion (1986)
#'
#' @param x a numeric matrix. Rows represent different variables, columns- observations.
#' @param mu_0 a numeric value. The scaling parameter of the target vector used in the
#' construction of the shrinkage estimator of the mean vector.
#' @return a numeric vector containing the Bayes-Stein shrinkage estimation of the mean vector
#' @references \insertRef{Jorion1986}{HDShOP}
#' @examples
#' n <- 7e2 # number of realizations
#' p <- .5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_bs(x=x, mu_0 = rep(1,p))
#' @export
mean_bs <- function(x, mu_0)
{
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)

  # Bayes-Stein mus and alphas
  alp_JS_hat <- as.numeric((p+2) / (p+2 + n*t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)))
  means <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect
  means
}


#' James-Stein shrinkage estimator of the mean vector
#'
#' James-Stein shrinkage estimator of the mean vector as suggested in Jorion (1986)
#'
#' @inheritParams mean_bs
#' @return a numeric vector containing the James-Stein shrinkage estimator of the mean vector
#' @references \insertRef{Jorion1986}{HDShOP}
#' @examples
#' n<-7e2 # number of realizations
#' p<-.5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_js(x=x, mu_0 = rep(1,p))
#' @export
mean_js <- function(x, mu_0)
{
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)

  # James-Stein mus and alphas
  val <- as.numeric( (p-2) / n / (t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)) )
  alp_JS_hat <- min(1,val)
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect
}

#' BOP shrinkage estimator
#'
#'Shrinkage estimator of the high-dimensional mean vector as suggested in Bodnar et al. (2019).
#'
#' @param x a numeric matrix. Rows represent different variables, columns- observations.
#' @param mu_0 a numeric vector. The target vector used in the construction of the shrinkage estimator of the mean vector.
#' @return a numeric vector containing the shrinkage estimation of the mean vector
#' @references \insertRef{BOP2019}{HDShOP}
#' @examples
#' n<-7e2 # number of realizations
#' p<-.5*n # number of assets
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#' mm <- mean_bop19(x=x, mu_0 = rep(1,p))
#' @export
mean_bop19 <- function(x, mu_0)
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
  output
}
