
#' constructor of EU portfolio object. IEEE 2020
#'
#' @inheritParams EUShrinkPortfolio
#' @param b a numeric value. The target for weight shrinkage.
#' @references \insertRef{BDOPS2020}{hdsp}
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_ExUtil_portfolio_weights_BDOPS20(x=x, gamma=gamma, b=b)
#' str(test)
#' @export
new_ExUtil_portfolio_weights_BDOPS20 <- function(x, gamma, b){

  p <- nrow(x)
  n <- ncol(x)

  if (is.data.frame(x)) x <- as.matrix(x)

  cov_mtrx <- Sigma_sample_estimator(x)
  means <- .rowMeans(x, m=p, n=n)

  invSS <- solve(cov_mtrx)
  Ip <- rep.int(1, nrow(x))
  Q_n_hat <- Q_hat_n(x) # this could be optimized

  W_EU_hat <- as.vector(
    (invSS %*% Ip)/as.numeric(t(Ip) %*% invSS %*% Ip) +
      Q_n_hat %*% means/gamma,
    mode = 'numeric')

  al <- alpha_hat_star_c(gamma, x=x, b=b)
  weights <- al*W_EU_hat + (1-al)*b
  structure(list(cov_mtrx=cov_mtrx,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 weights=weights,
                 alpha=al),
            class = c("ExUtil_portfolio_weights_BDOPS20", "ExUtil_portfolio")) # add alpha, stand dev, p-value when type=weights
}


