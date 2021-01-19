#### Direct cov


#' Ledoit and Wolf 2020
#'
#' @inheritParams EUShrinkPortfolio
#' @references \insertRef{LW2020}{hdsp}
#' @export
new_ExUtil_portfolio_cov_LW02 <- function(x, gamma){

  if (is.data.frame(x)) x <- as.matrix(x)

  cov_mtrx <- nonlin_shrinkLW(x)
  invSS <- solve(cov_mtrx)

  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  structure(list(cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio_cov_LW02", "ExUtil_portfolio"))
}

#' Bodnar, Gupta, Parolya 2014
#'
#' @inheritParams EUShrinkPortfolio
#' @param TM the target matrix for shrinkage of covariance
#' @references \insertRef{BGP2014}{hdsp}
#' @export
new_ExUtil_portfolio_cov_BGP14 <- function(x, gamma, TM){

  p <- nrow(x)
  n <- ncol(x)

  if (is.data.frame(x)) x <- as.matrix(x)

  SCM <- Sigma_sample_estimator(x)
  cov_mtrx <- CovShrinkBGP14(n=n, TM=TM, SCM=SCM)
  invSS <- solve(cov_mtrx)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  structure(list(cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio_cov_BGP14", "ExUtil_portfolio"))
}


#### Inverse cov shrinkage
#' Portfolio with inverse covariance shrinkage
#'
#' @inheritParams EUShrinkPortfolio
#' @param TM the target matrix for shrinkage of the inverse covariance matrix
#' @references \insertRef{BGP2016}{hdsp}
#' @export
new_ExUtil_portfolio_icov_BGP16 <- function(x, gamma, TM){

  p <- nrow(x)
  n <- ncol(x)

  if (is.data.frame(x)) x <- as.matrix(x)

  SCM <- Sigma_sample_estimator(x)
  iSCM=solve(SCM)
  invSS <- InvCovShrinkBGP16(n=n, p=p, TM,iSCM)

  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  structure(list(inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio_pm_BGP16", "ExUtil_portfolio"))
}





