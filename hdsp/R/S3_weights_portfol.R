#' constructor of EU portfolio object. IEEE 2020
#'
#' @inheritParams EUShrinkPortfolio
#' @param b a numeric variable. The target for weight shrinkage.
#' @param alph a numeric variable. The level of confidence for weight intervals
#' @references \insertRef{BDOPS2020}{hdsp}
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- new_ExUtil_portfolio_weights_BDOPS20(x=x, gamma=gamma, b=b, alph=0.05)
#' str(test)
#' @export
new_ExUtil_portfolio_weights_BDOPS20 <- function(x, gamma, b, alph){

  cl <- match.call()
  p <- nrow(x)
  n <- ncol(x)
  cc<- p/n
  if (is.data.frame(x)) x <- as.matrix(x)


  #### Direct / inverse covariance computation
  cov_mtrx <- Sigma_sample_estimator(x)
  Message <- 'No problems occured'
  invSS <- tryCatch(solve(cov_mtrx), error = function(e) e)
  if(!is.matrix(invSS)) {
    invSS <- MASS::ginv(cov_mtrx)
    invSS <- invSS[1:p, 1:p]
    Message <- 'solve() failed, falled back to generalized inverse'
  }

  mu_est <- .rowMeans(x, m=p, n=n)
  ones <- rep.int(1, p)
  tones<-t(ones)

  ########
  V_hat_c <- V_hat_c_fast(ones=ones, invSS=invSS, tones=tones, c=cc)  # V_hat_c is V.est
  Q_n_hat <- Q_hat_n_fast(invSS=invSS, Ip=ones, tIp=tones)  # Q_n_hat is Q.est
  s_hat_c <- as.numeric((1-cc)*t(mu_est)%*% Q_n_hat %*% mu_est-cc) # s.est
  R_hat_GMV<-(tones%*% invSS %*% mu_est)/as.numeric(tones%*%invSS%*%ones) # R.est, returns of EU portfolio

  ######## calculate shrinkage weights for EU or GMVP ##############
  V_hat_b <- V_b(Sigma=cov_mtrx, b=b)  # Vb.est
  R_hat_b <- R_b(mu=mu_est, b=b) # Rb.est


  W_EU_hat <- as.vector(
    (invSS %*% ones)/as.numeric(tones %*% invSS %*% ones) +
      Q_n_hat %*% mu_est/gamma,
    mode = 'numeric')

  # alpha_EU
  al <- alpha_hat_star_c_fast(gamma=gamma, c=cc, s=s_hat_c, b=b, R_GMV=R_hat_GMV,
                              R_b=R_hat_b, V_c=V_hat_c, V_b=V_hat_b)
  weights <- al*W_EU_hat + (1-al)*b # w_EU_shr


  #### Confidence intervals for weights ####

  T_dens <- matrix(data=NA, nrow=p, ncol=5)

  for(i in seq(1,p,by=1)) {

    L <- matrix(rep(0,p), nrow=1)
    L[1,i]<- 1 # select one!!! component of weights vector for a test
    eta.est<-((s_hat_c+cc)/s_hat_c)*L%*%Q_n_hat%*%mu_est/as.numeric(t(mu_est)%*%Q_n_hat%*%mu_est)
    w.est<-(L%*%invSS%*%ones)/sum(tones%*%invSS%*%ones)+gamma^{-1}*s_hat_c*eta.est # weight (component of weights vector)

    ############# BDOPS, formula 15
    Omega.Lest<- Omega.Lest(s_hat_c=s_hat_c, cc=cc, gamma=gamma,
                            V_hat_c=V_hat_c, L=L, Q_n_hat=Q_n_hat, eta.est=eta.est)

    TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est)
    low_bound <- w.est - qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
    upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
    p_value<- pchisq(TDn, df=1, lower.tail=FALSE)

    T_dens[i,] <- c(weights[i], low_bound, upp_bound, TDn,  p_value) #first column= shrinkage estimator for EU Portfolio weights
  }

  colnames(T_dens) <- c('weight', 'low_bound', 'upp_bound', 'TDn',  'p_value')

  Port_Var <- as.numeric(t(weights)%*%cov_mtrx%*%weights)
  Port_mean_return <- as.numeric(mu_est %*% weights)
  Sharpe <- Port_mean_return/sqrt(Port_Var)

  structure(list(call=cl,
                 cov_mtrx=cov_mtrx,
                 inv_cov_mtrx=invSS,
                 means=mu_est,
                 W_EU_hat=W_EU_hat,
                 weights=weights,
                 alpha=al,
                 Port_Var=Port_Var,
                 Port_mean_return=Port_mean_return,
                 Sharpe=Sharpe,
                 weight_intervals=T_dens,
                 Message=Message),
            class = c("ExUtil_portfolio_weights_BDOPS20", "ExUtil_portfolio"))
  }


#' constructor of GMVA portfolio object.
#'
#' @inheritParams EUShrinkPortfolio
#' @param b a numeric value. The target for weight shrinkage.
#' @inheritParams new_ExUtil_portfolio_weights_BDOPS20
#' @references \insertRef{BDOPS2020}{hdsp}
#' @examples
#' library(MASS)
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#'
#' Mtrx <- RandCovMtrx(n=n, p=p, q=20.55, mu=seq(0.2,-0.2, length.out=p))
#' x <- t(mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))
#'
#' test <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, alph=0.05)
#' str(test)
#'
#' @export
new_GMV_portfolio_weights_BDPS19 <- function(x, b, alph){

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
  Message <- 'No problems occured'
  iS <- tryCatch(solve(cov_mtrx), error = function(e) e)
  if(!is.matrix(iS)) {
    iS <- MASS::ginv(cov_mtrx)
    iS <- iS[1:p, 1:p]
    Message <- 'solve() failed, falled back to generalized inverse'
  }

  V.est <- V_hat_c_fast(ones=ones, invSS=iS, tones=tones, c=cc)
  # V.est<-(1-cc)^{-1}/sum(tones%*%iS%*%ones) #estimated variance (cons.)
  # Q.est<- iS-(iS%*%ones%*%tones%*%iS)/sum(tones%*%iS%*%ones)
  Q.est <- Q_hat_n_fast(invSS=iS, Ip=ones, tIp=tones)

  ####  for calculating shrinkage GMVP weights
  w_GMVP_whole <- iS%*%ones/as.numeric(tones%*%iS%*%ones) # sample estimator for GMVP weights
  Lb <- (1-cc)*t(b)%*%cov_mtrx%*%b*as.numeric(tones%*%iS%*%ones) -1 #relative loss of GMVP f.(16) IEEE 2019
  alpha_GMVP <- as.numeric((1-cc)*Lb/(cc+(1-cc)*Lb)) #shrinkage intensity f.(16) IEEE 2019
  w_GMV_shr <- alpha_GMVP*w_GMVP_whole + (1-alpha_GMVP)*b # f.(17) IEEE 2019


  #### Confidence intervals for weights ####

  T_dens <- matrix(data=NA, nrow=p, ncol=5)

  for(i in seq(1,p,by=1)) {

    L <- matrix(rep(0,p), nrow=1)
    L[1,i]<- 1 # select one!!! component of weights vector for a test
    w.est <- (L%*%iS%*%ones)/sum(tones%*%iS%*%ones) #(16) for GMVP
    ############# BDOPS, formula 15
    Omega.Lest <- V.est*(1-cc)*L%*%Q.est%*%t(L) # (17) IEEE 2021 simplifies to that)

    TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est) #test statistics

    low_bound <- w.est - qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
    upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)

    p_value <- pchisq(TDn, df=1, lower.tail=FALSE)

    #return(c(w.est, low_bound, upp_bound, TDn,  p_value)) #if you want to have the components of the weights vector in F. ()
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
                 weight_intervals=T_dens,
                 Message=Message),
            class = c("GMV_portfolio_weights_BDPS19", "ExUtil_portfolio"))
}
