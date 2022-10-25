
#' Mean-variance shrinkage portfolio when c > 1.
#'
#' TBA
#'
#' @examples
#' # Assets with a diagonal covariance matrix
#'
#' n<-3e2 # number of realizations
#' p<-1.3*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_MV_portfolio_weights_BDOPS21_pgn(x=x, gamma=gamma, b=b, beta=0.05)
#' summary(test)
#'
#' # Assets with a non-diagonal covariance matrix
#'
#' Mtrx <- RandCovMtrx(p=p)
#' x <- t(MASS::mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))
#'
#' test <- new_MV_portfolio_weights_BDOPS21_pgn(x=x, gamma=gamma, b=b, beta=0.05)
#' str(test)
#'
#' @export
new_MV_portfolio_weights_BDOPS21_pgn <- function(x, gamma, b, beta){

  cl<- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc<- p/n
  if (is.data.frame(x)) x <- as.matrix(x)


  #### Direct / inverse covariance computation
  cov_mtrx <- HDShOP:::Sigma_sample_estimator(x)
  invSS <- MASS::ginv(cov_mtrx)

  mu_est <- .rowMeans(x, m=p, n=n)
  ones <- rep.int(1, p)
  tones <- t(ones)

  ########
  V_hat_c_pgn <- HDShOP:::V_hat_c_pgn_fast(ones=ones, invSS=invSS, tones=tones, c=cc)  #
  Q_n_pgn_hat <- HDShOP:::Q_hat_n_fast(invSS=invSS, Ip=ones, tIp=tones)  # Q_n_hat is Q.est
  s_hat_c_pgn <- as.numeric(cc*(cc-1)*t(mu_est) %*% Q_n_pgn_hat %*% mu_est-cc) # s.est
  R_hat_GMV   <- (tones%*% invSS %*% mu_est)/as.numeric(tones%*%invSS%*%ones) # R.est, returns of EU portfolio

  ######## calculate shrinkage weights for EU or GMVP ##############
  V_hat_b <- HDShOP:::V_b(Sigma=cov_mtrx, b=b)  # Vb.est
  R_hat_b <- HDShOP:::R_b(mu=mu_est, b=b) # Rb.est


  W_EU_hat <- as.vector(
      (invSS %*% ones)/as.numeric(tones %*% invSS %*% ones) +
      Q_n_pgn_hat %*% mu_est/gamma,
    mode = 'numeric')

  # alpha_EU
  al <- HDShOP:::alpha_hat_star_c_pgn_fast(gamma=gamma, c=cc, s=s_hat_c_pgn, R_GMV=R_hat_GMV,
                                  R_b=R_hat_b, V_c=V_hat_c_pgn, V_b=V_hat_b)
  weights <- al*W_EU_hat + (1-al)*b # w_EU_shr


  #### Confidence intervals for weights are omitted ####

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
                 Sharpe=Sharpe),
            class = c("MeanVar_portfolio"))
}




#' Global minimum-variance shrinkage portfolio when c > 1.
#'
#' TBA
#'
#' @examples
#'
#' n<-3e2 # number of realizations
#' p<-1.3*n # number of assets
#' b<-rep(1/p,p)
#'
#' # Assets with a diagonal covariance matrix
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_GMV_portfolio_weights_BDPS19_pgn(x=x, b=b, beta=0.05)
#' str(test)
#'
#' # Assets with a non-diagonal covariance matrix
#' Mtrx <- RandCovMtrx(p=p)
#' x <- t(MASS::mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))
#'
#' test <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, beta=0.05)
#' summary(test)
#'
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
  cov_mtrx <- HDShOP:::Sigma_sample_estimator(x)
  iS <- MASS::ginv(cov_mtrx)

  ########
  V_hat_c_pgn <- HDShOP:::V_hat_c_pgn_fast(ones=ones, invSS=iS, tones=tones, c=cc)  #

  ####  for calculating shrinkage GMVP weights
  V_hat_b <- HDShOP:::V_b(Sigma=cov_mtrx, b=b)  # Vb.est

  w_GMVP_whole <- as.vector(iS%*%ones/as.numeric(tones%*%iS%*%ones),
                            mode='numeric') # sample estimator for GMVP weights

  alpha_GMVP <- HDShOP:::alpha_hat_star_c_GMV_pgn_fast(c=cc, V_c=V_hat_c_pgn, V_b=V_hat_b)
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










