# Direct cov


# Ledoit and Wolf 2002
new_ExUtil_portfolio_cov_LW02 <- function(x, gamma){

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
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio","ExUtil_portfolio_cov_LW02"))
}

# Bodnar, Gupta, Parolya 2014
new_ExUtil_portfolio_cov_BGP14 <- function(x, gamma, TM){

  SCM <- Sigma_sample_estimator(x)
  cov_mtrx <- CovShrinkBGP(TM=TM, SCM=SCM)
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
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio","ExUtil_portfolio_cov_BGP14"))
}


# Inverse cov shrinkage

new_ExUtil_portfolio_pm_BGP16 <- function(x, gamma, TM){

  SCM <- Sigma_sample_estimator(x)
  iSCM=solve(SCM)
  invSS <- InvCovShrinkBGP(TM,iSCM)

  p <- nrow(x)
  n <- ncol(x)
  means <- .rowMeans(x, m=p, n=n)
  I_vect <- rep(1, times=p)
  Q_n_hat <- invSS - (invSS %*% I_vect %*% t(I_vect) %*% invSS)/as.numeric(t(I_vect) %*% invSS %*% I_vect)

  W_EU_hat <- as.vector(
    (invSS %*% I_vect)/as.numeric(t(I_vect) %*% invSS %*% I_vect) +
      Q_n_hat %*% means / gamma,
    mode = 'numeric')

  structure(list(cov_mtrx=solve(invSS),
                 invSS=invSS,
                 means=means,
                 W_EU_hat=W_EU_hat),
            class = c("ExUtil_portfolio","ExUtil_portfolio_pm_BGP16"))
}





