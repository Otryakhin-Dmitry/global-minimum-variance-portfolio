#' @export
new_ExUtil_portfolio_weights_BDOPS20_new <- function(x, gamma, b, alph){

  p <- nrow(x)
  n <- ncol(x)
  cc<- p/n

  if (is.data.frame(x)) x <- as.matrix(x)

  cov_mtrx <- Sigma_sample_estimator(x)
  mu_est <- .rowMeans(x, m=p, n=n)

  invSS <- solve(cov_mtrx)
  ones <- rep.int(1, nrow(x))
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

  structure(list(cov_mtrx=cov_mtrx,
                 means=mu_est,
                 W_EU_hat=W_EU_hat,
                 weights=weights,
                 alpha=al,
                 weight_intervals=T_dens),
            class = c("ExUtil_portfolio_weights_BDOPS20", "ExUtil_portfolio"))
}
