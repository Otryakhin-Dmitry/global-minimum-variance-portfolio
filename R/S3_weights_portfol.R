#' Constructor of MV portfolio object
#'
#' Constructor of mean-variance shrinkage portfolios.
#' new_MV_portfolio_weights_BDOPS21 is for the case p<n, while
#' new_MV_portfolio_weights_BDOPS21_pgn is for p>n, where p is the number of
#' assets and n is the number of observations.
#' For more details of the method, see \code{\link{MVShrinkPortfolio}}.
#'
#' @inheritParams MVShrinkPortfolio
#' @param b a numeric variable. The weights of the target portfolio.
#' @param beta a numeric variable. The confidence level for weight intervals.
#' @return an object of class MeanVar_portfolio with subclass
#' MV_portfolio_weights_BDOPS21.
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the asset returns |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean vector of the asset returns |
#' | W_mv_hat | sample estimator of the portfolio weights |
#' | weights | shrinkage estimator of the portfolio weights |
#' | alpha | shrinkage intensity for the weights |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | expected portfolio return |
#' | Sharpe | portfolio Sharpe ratio |
#' | weight_intervals | A data frame, see details |
#'
#' weight_intervals contains a shrinkage estimator of portfolio weights,
#' asymptotic confidence intervals for the true portfolio weights, value of
#' the test statistic and the p-value of the test on the equality of
#' the weight of each individual asset to zero
#' (see Section 4.3 of Bodnar, Dette, Parolya and Thorsén 2021).
#' weight_intervals is only computed when p<n.
#' @md
#'
#' @references \insertRef{BDOPS2021}{HDShOP}
#' @references \insertRef{BDNT21}{HDShOP}
#' @examples
#'
#' # c<1
#'
#' # Assets with a diagonal covariance matrix
#'
#' n <- 3e2 # number of realizations
#' p <- .5*n # number of assets
#' b <- rep(1/p,p)
#' gamma <- 1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=gamma, b=b, beta=0.05)
#' summary(test)
#'
#' # Assets with a non-diagonal covariance matrix
#'
#' Mtrx <- RandCovMtrx(p=p)
#' x <- t(MASS::mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))
#'
#' test <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=gamma, b=b, beta=0.05)
#' str(test)
#'
#'
#' # c>1
#'
#' n <-2e2 # number of realizations
#' p <-1.2*n # number of assets
#' b <-rep(1/p,p)
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_MV_portfolio_weights_BDOPS21_pgn(x=x, gamma=gamma,
#'                                              b=b, beta=0.05)
#' summary(test)
#'
#' # Assets with a non-diagonal covariance matrix
#'
#' @export
new_MV_portfolio_weights_BDOPS21 <- function(x, gamma, b, beta){

  cl <- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc <- p/n
  if (is.data.frame(x)) x <- as.matrix(x)


  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- solve(cov_mtrx)

  mu_est <- .rowMeans(x, m=p, n=n)
  ones <- rep.int(1, p)
  tones <- t(ones)

  ########
  # V_hat_c is V.est
  V_hat_c <- V_hat_c_fast(ones=ones, invSS=invSS, tones=tones, c=cc)
  Q_n_hat <- Q_hat_n_fast(invSS=invSS, Ip=ones, tIp=tones)  # Q_n_hat is Q.est
  s_hat_c <- as.numeric((1-cc)*t(mu_est)%*% Q_n_hat %*% mu_est-cc) # s.est
  # R.est, returns of EU portfolio
  R_hat_GMV <- (tones%*% invSS %*% mu_est)/as.numeric(tones%*%invSS%*%ones)

  ######## calculate shrinkage weights for EU or GMVP ##############
  V_hat_b <- V_b(Sigma=cov_mtrx, b=b)  # Vb.est
  R_hat_b <- R_b(mu=mu_est, b=b) # Rb.est


  W_EU_hat <- as.vector(
    (invSS %*% ones)/as.numeric(tones %*% invSS %*% ones) +
      Q_n_hat %*% mu_est/gamma,
    mode = 'numeric')

  # alpha_EU
  al <- alpha_hat_star_c_fast(gamma=gamma, c=cc, s=s_hat_c, R_GMV=R_hat_GMV,
                              R_b=R_hat_b, V_c=V_hat_c, V_b=V_hat_b)
  weights <- al*W_EU_hat + (1-al)*b # w_EU_shr


  #### Confidence intervals for weights ####

  T_dens <- matrix(data=NA, nrow=p, ncol=5)

  for(i in seq(1,p,by=1)) {

    L <- matrix(rep(0,p), nrow=1)
    L[1,i] <- 1 # select one!!! component of weights vector for a test
    eta.est <- ((s_hat_c+cc)/s_hat_c)*L%*%Q_n_hat%*%mu_est/
               as.numeric(t(mu_est)%*%Q_n_hat%*%mu_est)
    w.est <- (L%*%invSS%*%ones)/sum(tones%*%invSS%*%ones) +
             gamma^{-1}*s_hat_c*eta.est # weight (component of weights vector)

    ############# BDOPS, formula 15
    Omega.Lest <- Omega.Lest(s_hat_c=s_hat_c, cc=cc, gamma=gamma,
                             V_hat_c=V_hat_c, L=L, Q_n_hat=Q_n_hat,
                             eta.est=eta.est)

    TDn <- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est)
    low_bound <- w.est - qnorm(1-beta/2)*sqrt(Omega.Lest)/sqrt(n-p)
    upp_bound <- w.est + qnorm(1-beta/2)*sqrt(Omega.Lest)/sqrt(n-p)
    p_value <- pchisq(TDn, df=1, lower.tail=FALSE)
    #first column= shrinkage estimator for EU Portfolio weights
    T_dens[i,] <- c(weights[i], low_bound, upp_bound, TDn,  p_value)
  }

  colnames(T_dens) <- c('weight', 'low_bound', 'upp_bound', 'TDn',  'p_value')

  Port_Var <- as.numeric(t(weights)%*%cov_mtrx%*%weights)
  Port_mean_return <- as.numeric(mu_est %*% weights)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=mu_est,
                 W_mv_hat=W_EU_hat,
                 weights=weights,
                 alpha=al,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe,
                 weight_intervals=T_dens),
            class = c("MV_portfolio_weights_BDOPS21", "MeanVar_portfolio"))
  }


#' @rdname new_MV_portfolio_weights_BDOPS21
#' @export
new_MV_portfolio_weights_BDOPS21_pgn <- function(x, gamma, b, beta){

  cl <- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc <- p/n
  if (is.data.frame(x)) x <- as.matrix(x)


  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  invSS <- MASS::ginv(cov_mtrx)

  mu_est <- .rowMeans(x, m=p, n=n)
  ones <- rep.int(1, p)
  tones <- t(ones)

  ########
  V_hat_c_pgn <- V_hat_c_pgn_fast(ones=ones, invSS=invSS, tones=tones, c=cc)
  # Q_n_hat is Q.est
  Q_n_pgn_hat <- Q_hat_n_fast(invSS=invSS, Ip=ones, tIp=tones)
  # s.est
  s_hat_c_pgn <- as.numeric(cc*(cc-1)*t(mu_est) %*% Q_n_pgn_hat %*% mu_est-cc)
  # R.est, returns of EU portfolio
  R_hat_GMV   <- (tones%*% invSS %*% mu_est)/as.numeric(tones%*%invSS%*%ones)

  ######## calculate shrinkage weights for EU or GMVP ##############
  V_hat_b <- V_b(Sigma=cov_mtrx, b=b)  # Vb.est
  R_hat_b <- R_b(mu=mu_est, b=b) # Rb.est


  W_EU_hat <- as.vector(
    (invSS %*% ones)/as.numeric(tones %*% invSS %*% ones) +
      Q_n_pgn_hat %*% mu_est/gamma,
    mode = 'numeric')

  # alpha_EU
  al <- alpha_hat_star_c_pgn_fast(gamma=gamma, c=cc, s=s_hat_c_pgn,
                                  R_GMV=R_hat_GMV, R_b=R_hat_b,
                                  V_c=V_hat_c_pgn, V_b=V_hat_b)
  weights <- al*W_EU_hat + (1-al)*b # w_EU_shr


  #### Confidence intervals for weights are omitted ####

  Port_Var <- as.numeric(t(weights)%*%cov_mtrx%*%weights)
  Port_mean_return <- as.numeric(mu_est %*% weights)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call = cl,
                 cov_mtrx = cov_mtrx,
                 inv_cov_mtrx = invSS,
                 means = mu_est,
                 W_mv_hat = W_EU_hat,
                 weights = weights,
                 alpha = al,
                 Port_Var = Port_Var,
                 Port_mean_return = Port_mean_return,
                 Sharpe = Sharpe),
            class = c("MeanVar_portfolio"))
}


#' Constructor of GMV portfolio object.
#'
#' Constructor of global minimum variance portfolio.
#' new_GMV_portfolio_weights_BDPS19 is for the case p<n, while
#' new_GMV_portfolio_weights_BDPS19_pgn is for p>n, where p is the number of
#' assets and n is the number of observations. For more details
#' of the method, see \code{\link{MVShrinkPortfolio}}.
#'
#' @inheritParams MVShrinkPortfolio
#' @param b a numeric vector. The weights of the target portfolio.
#' @inheritParams new_MV_portfolio_weights_BDOPS21
#' @return an object of class MeanVar_portfolio with subclass
#' GMV_portfolio_weights_BDPS19.
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the asset returns |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean vector estimate of the asset returns |
#' | w_GMVP | sample estimator of portfolio weights |
#' | weights | shrinkage estimator of portfolio weights |
#' | alpha | shrinkage intensity for the weights |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | expected portfolio return |
#' | Sharpe | portfolio Sharpe ratio |
#' | weight_intervals | A data frame, see details |
#'
#' weight_intervals contains a shrinkage estimator of portfolio weights,
#' asymptotic confidence intervals for the true portfolio weights,
#' the value of test statistic and the p-value of the test on the equality of
#' the weight of each individual asset to zero
#' \insertCite{@see Section 4.3 of @BDNT21}{HDShOP}.
#' weight_intervals is only computed when p<n.
#' @md
#'
#' @references \insertRef{BDPS2019}{HDShOP}
#' @references \insertRef{BPS2018}{HDShOP}
#' @references \insertRef{BDNT21}{HDShOP}
#' @examples
#'
#' # c<1
#'
#' n <- 3e2 # number of realizations
#' p <- .5*n # number of assets
#' b <- rep(1/p,p)
#'
#' # Assets with a diagonal covariance matrix
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, beta=0.05)
#' str(test)
#'
#' # Assets with a non-diagonal covariance matrix
#' Mtrx <- RandCovMtrx(p=p)
#' x <- t(MASS::mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))
#'
#' test <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, beta=0.05)
#' summary(test)
#'
#' # c>1
#'
#' p <- 1.3*n # number of assets
#' b <- rep(1/p,p)
#'
#' # Assets with a diagonal covariance matrix
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_GMV_portfolio_weights_BDPS19_pgn(x=x, b=b, beta=0.05)
#' str(test)
#'
#' @export
new_GMV_portfolio_weights_BDPS19 <- function(x, b, beta){

  cl <- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc <- p/n
  if (is.data.frame(x)) x <- as.matrix(x)

  ones <- rep.int(1, p)
  tones <- t(ones)
  mu_est <- .rowMeans(x, m=p, n=n)


  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  iS <- solve(cov_mtrx)

  V.est <- V_hat_c_fast(ones=ones, invSS=iS, tones=tones, c=cc)
  Q.est <- Q_hat_n_fast(invSS=iS, Ip=ones, tIp=tones)

  ####  for calculating shrinkage GMVP weights
  w_GMVP_whole <- as.vector(iS%*%ones/as.numeric(tones%*%iS%*%ones),
                            mode='numeric') # sample estimator for GMVP weights
  #relative loss of GMVP f.(16) IEEE 2019
  Lb <- (1-cc)*t(b)%*%cov_mtrx%*%b*as.numeric(tones%*%iS%*%ones) -1
  #shrinkage intensity f.(16) IEEE 2019
  alpha_GMVP <- as.numeric((1-cc)*Lb/(cc+(1-cc)*Lb))
  w_GMV_shr <- alpha_GMVP*w_GMVP_whole + (1-alpha_GMVP)*b # f.(17) IEEE 2019


  #### Confidence intervals for weights ####

  T_dens <- matrix(data=NA, nrow=p, ncol=5)

  for(i in seq(1,p,by=1)) {

    L <- matrix(rep(0,p), nrow=1)
    L[1,i]<- 1 # select one!!! component of weights vector for a test
    w.est <- (L%*%iS%*%ones)/sum(tones%*%iS%*%ones) #(16) for GMVP
    ############# BDOPS, formula 15
    # (17) IEEE 2021 simplifies to that)
    Omega.Lest <- V.est*(1-cc)*L%*%Q.est%*%t(L)

    TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est) #test statistics

    low_bound <- w.est - qnorm(1-beta/2)*sqrt(Omega.Lest)/sqrt(n-p)
    upp_bound <- w.est + qnorm(1-beta/2)*sqrt(Omega.Lest)/sqrt(n-p)

    p_value <- pchisq(TDn, df=1, lower.tail=FALSE)

    T_dens[i,] <- c(w_GMV_shr[i], low_bound, upp_bound, TDn,  p_value)
  }

  colnames(T_dens) <- c('weight', 'low_bound', 'upp_bound', 'TDn',  'p_value')

  Port_Var <- 1/as.numeric(tones%*%iS%*%ones)
  Port_mean_return <- as.numeric(mu_est %*% w_GMV_shr)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  #### Output
  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=iS,
                 means=mu_est,
                 w_GMVP=w_GMVP_whole,
                 weights=w_GMV_shr,
                 alpha=alpha_GMVP,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe,
                 weight_intervals=T_dens),
            class = c("GMV_portfolio_weights_BDPS19", "MeanVar_portfolio"))
}


#' @rdname new_GMV_portfolio_weights_BDPS19
#' @export
new_GMV_portfolio_weights_BDPS19_pgn <- function(x, b, beta){

  cl <- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc<- p/n
  if (is.data.frame(x)) x <- as.matrix(x)

  ones <- rep.int(1, p)
  tones<-t(ones)
  mu_est <- .rowMeans(x, m=p, n=n)


  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  iS <- MASS::ginv(cov_mtrx)

  ########
  V_hat_c_pgn <- V_hat_c_pgn_fast(ones=ones, invSS=iS, tones=tones, c=cc)  #

  ####  for calculating shrinkage GMVP weights
  V_hat_b <- V_b(Sigma=cov_mtrx, b=b)  # Vb.est

  w_GMVP_whole <- as.vector(iS%*%ones/as.numeric(tones%*%iS%*%ones),
                            mode='numeric') # sample estimator for GMVP weights

  alpha_GMVP <- alpha_hat_star_c_GMV_pgn_fast(c=cc, V_c=V_hat_c_pgn,
                                              V_b=V_hat_b)
  w_GMV_shr <- alpha_GMVP*w_GMVP_whole + (1-alpha_GMVP)*b


  #### Confidence intervals for weights are omitted ####

  Port_Var <- 1/as.numeric(tones%*%iS%*%ones)
  Port_mean_return <- as.numeric(mu_est %*% w_GMV_shr)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  #### Output
  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=iS,
                 means=mu_est,
                 w_GMVP=w_GMVP_whole,
                 weights=w_GMV_shr,
                 alpha=alpha_GMVP,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe),
            class = c("MeanVar_portfolio"))
}
