
#' constructor of EU portfolio object. Type=mean, Mean.type=Bayes-Stein.
#'
#' @inheritParams EUShrinkPortfolio
#' @param mu_0 a numeric value. The scaling of the target for shrinkage of the mean vector.
#' @references \insertRef{Jorion1986}{hdsp}
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

  if (is.data.frame(x)) x <- as.matrix(x)

  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  # Bayes-Stein mus and alphas
  alp_JS_hat <- as.numeric((p+2) / (p+2 + n*t(means-mu_0*I_vect)%*%invSS%*%(means-mu_0*I_vect)))
  mu_hat_JS <- (1-alp_JS_hat) * means + alp_JS_hat * mu_0 * I_vect

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
    Q_n_hat %*% mu_hat_JS / gamma,
    mode = 'numeric')

  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 mu_hat_JS=mu_hat_JS,
                 alp_JS_hat=alp_JS_hat,
                 W_EU_hat=W_EU_hat
                 ),
            class = c("ExUtil_portfolio","ExUtil_portfolio_mean_Bayes-Stein")) # add alpha, stand dev, p-value when type=weights
}


#' constructor of EU portfolio object. Type=mean, Mean.type=James-Stein.
#'
#' @inheritParams new_ExUtil_portfolio_mean_BayesStein
#' @references \insertRef{Jorion1986}{hdsp}
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

  if (is.data.frame(x)) x <- as.matrix(x)

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

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
    Q_n_hat %*% mu_hat_JS / gamma,
    mode = 'numeric')

  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 mu_hat_JS=mu_hat_JS,
                 alp_JS_hat=alp_JS_hat,
                 W_EU_hat=W_EU_hat
                 ),
            class = c("ExUtil_portfolio","ExUtil_portfolio_mean_James-Stein")) # add alpha, stand dev, p-value when type=weights
}


#' constructor of EU portfolio object. Type=mean, Mean.type=BOP.
#'
#' @inheritParams EUShrinkPortfolio
#' @param mu_0 a numeric vector. The target for shrinkage of the mean vector of asset returns.
#' @references \insertRef{BOP2019}{hdsp}
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

  if (is.data.frame(x)) x <- as.matrix(x)

  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)
  p <- nrow(x)
  n <- ncol(x)
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

  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 mu_hat_BOP=mu_hat_BOP,
                 alpha_star_hat=alpha_star_hat,
                 beta_star_hat=beta_star_hat,
                 W_EU_hat=W_EU_hat
  ),
  class = c("ExUtil_portfolio","ExUtil_portfolio_mean_BOP")) # add alpha, stand dev, p-value when type=weights
}
