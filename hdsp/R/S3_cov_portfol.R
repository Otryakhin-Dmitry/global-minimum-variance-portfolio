#### Direct cov


#' constructor of EU portfolio object. Ledoit and Wolf 2020
#'
#' @inheritParams EUShrinkPortfolio
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the covariance matrix estimate of the assets according to Ledoit and Wolf 2020|
#' | inv_cov_mtrx | the inverse of the covariance matrix estimate |
#' | means | sample mean vector estimate for the assets |
#' | W_EU_hat | portfolio weights estimate obtained via the above cov_mtrx and means |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @references \insertRef{LW2020}{hdsp}
#' @examples
#' c <- 0.4
#' n <- 7e2
#' p <- c*n
#' gamma <- 1
#'
#' X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)
#' port_LW02 <- new_ExUtil_portfolio_cov_LW02(x=X, gamma)
#' str(port_LW02)
#'
#' @export
new_ExUtil_portfolio_cov_LW02 <- function(x, gamma){

  cl <- match.call()
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

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- means %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("ExUtil_portfolio_cov_LW02", "ExUtil_portfolio"))
}

#' constructor of EU portfolio object. Bodnar, Gupta, Parolya 2014
#'
#' @inheritParams EUShrinkPortfolio
#' @param TM the target matrix for shrinkage of covariance
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the covariance matrix estimate of the assets according to Bodnar et al, 2014|
#' | inv_cov_mtrx | the inverse of the covariance matrix estimate |
#' | means | sample mean vector estimate for the assets |
#' | W_EU_hat | portfolio weights estimate obtained via the above cov_mtrx and means |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @references \insertRef{BGP2014}{hdsp}
#' @examples
#' c <- 0.4
#' n <- 7e2
#' p <- c*n
#' gamma <- 1
#'
#' Sigma=matrix(0, p, p)
#' diag(Sigma) <- 1
#'
#' X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)
#' port_BGP14 <- new_ExUtil_portfolio_cov_BGP14(x=X, gamma, TM=Sigma)
#' str(port_BGP14)
#'
#' @export
new_ExUtil_portfolio_cov_BGP14 <- function(x, gamma, TM){

  cl <- match.call()
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

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%cov_mtrx%*%W_EU_hat
  Port_mean_return <- means %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("ExUtil_portfolio_cov_BGP14", "ExUtil_portfolio"))
}


#### Inverse cov shrinkage
#' Portfolio with inverse covariance shrinkage
#'
#' @inheritParams EUShrinkPortfolio
#' @param TM the target matrix for shrinkage of the inverse covariance matrix
#' @return an object of class ExUtil_portfolio with a subclass ... .
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the assets |
#' | inv_cov_mtrx | the inverse of the covariance matrix estimate obtained via the method from Bodnar et al, 2016 |
#' | means | sample mean vector estimate for the assets |
#' | W_EU_hat | portfolio weights estimate obtained via the above inv_cov_mtrx and means |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | portfolio mean returns |
#' | Sharpe | portfolio Sharpe coefficient |
#' @md
#'
#' @references \insertRef{BGP2016}{hdsp}
#' @examples
#' c <- 0.4
#' n <- 7e2
#' p <- c*n
#' gamma <- 1
#'
#' Sigma=matrix(0, p, p)
#' diag(Sigma) <- 1
#'
#' X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)
#' port_BGP16 <- new_ExUtil_portfolio_icov_BGP16(x=X, gamma, TM=Sigma)
#' str(port_BGP16)
#'
#' @export
new_ExUtil_portfolio_icov_BGP16 <- function(x, gamma, TM){

  cl <- match.call()
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

  # Sharpe, mean return and portfolio variance
  Port_Var <- t(W_EU_hat)%*%SCM%*%W_EU_hat
  Port_mean_return <- means %*% W_EU_hat
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=SCM,
                 inv_cov_mtrx=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("ExUtil_portfolio_pm_BGP16", "ExUtil_portfolio"))
}





