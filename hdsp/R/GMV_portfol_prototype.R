
#' @export
new_GMV_portfolio_weights_BDPS19 <- function(x, gamma, b, alph){

  p <- ncol(x)
  n <- nrow(x)
  cc<- p/n
  #b <-rep(1/p, p) #for shrinkage, target

  ones<-matrix(1,p,1)
  tones<-t(ones)


  ####
  S <- cov(x)
  mu_est<-colMeans(x)
  iS<-solve(S)
  V.est<-(1-cc)^{-1}/sum(tones%*%iS%*%ones) #estimated variance (cons.)
  Q.est<- iS-(iS%*%ones%*%tones%*%iS)/sum(tones%*%iS%*%ones)


  ####  for calculating shrinkage GMVP weights
  w_GMVP_whole<- iS%*%ones/sum(tones%*%iS%*%ones) # sample estimator for GMVP weights
  Lb<- (1-cc)*t(b)%*%S%*%b*sum(tones%*%iS%*%ones) -1 #relative loss of GMVP f.(16) IEEE 2019
  alpha_GMVP<- as.numeric((1-cc)*Lb/(cc+(1-cc)*Lb)) #shrinkage intensity f.(16) IEEE 2019
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


  #### Output
  structure(list(cov_mtrx=S,
                 means=mu_est,
                 w_GMVP=w_GMVP_whole,
                 weights=w_GMV_shr,
                 alpha=alpha_GMVP,
                 weight_intervals=T_dens),
            class = c("GMV_portfolio_weights_BDPS19", "ExUtil_portfolio"))
}
