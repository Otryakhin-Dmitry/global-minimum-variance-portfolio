
library("doParallel")
no_cores <- detectCores()
library("foreach")
cl<-makeCluster(no_cores)
registerDoParallel(cl)



test_each_asset <- function(gamma, X_asset, alph){

  p <- ncol(X_asset)
  n <- nrow(X_asset)
  cc<- p/n
  b <-rep(1/p, p) #for shrinkage, target

  ones<-matrix(1,p,1)
  tones<-t(ones)

  S <- cov(X_asset)
  mu_est<-colMeans(X_asset)
  iS<-solve(S)
  V.est<-(1-cc)^{-1}/sum(tones%*%iS%*%ones) #estimated variance (cons.)
  Q.est<- iS-(iS%*%ones%*%tones%*%iS)/sum(tones%*%iS%*%ones)
  s.est<- as.numeric((1-cc)*t(mu_est)%*% Q.est %*% mu_est-cc)
  R.est<- (tones%*% iS %*% mu_est)/sum(tones%*%iS%*%ones) #returns of EU portfolio

  ######## calculate shrinkage weights for EU or GMVP ##############
  Vb.est<-t(b)%*%S%*%b     #optional:can be removed if you decide to import an existing function for shrinakge portolio weights
  Rb.est<-t(b) %*% mu_est

  #########################

  #########################
  ####  for calculating shrinkage EU weights (can be remove if you import your own function)
  w_EU_whole<- iS%*%ones/sum(tones%*%iS%*%ones) +
               gamma^{-1}*Q.est%*%mu_est # sample estimator for EU weights f.(6) IEEE 2021

  alpha_EU <- alpha_hat_star_c_fast(gamma=gamma, c=cc, s=s.est, b=b, R_GMV, R_b, V_GMV, V_b)
  alpha_EU <- as.numeric(gamma^{-1}*((R.est-Rb.est)*(1+1/(1-cc)) + gamma*(Vb.est-V.est)+gamma^{-1}*s.est/(1-cc))/(V.est/(1-cc)-2*(V.est+gamma^{-1}*(Rb.est-R.est)/(1-cc))+ gamma^{-2}*(s.est/(1-cc)^3+cc/(1-cc)^3)+Vb.est))
  #f.(25) IEEE 21

  w_EU_shr <- alpha_EU*w_EU_whole + (1-alpha_EU)*b # f.(26) IEEE 2021
  #end
  #######################


  T.dens = foreach(i=1:p, .combine="rbind", .inorder=FALSE, .packages=c("MASS")) %dopar%
    {
      ## EU case
        L <- matrix(rep(0,p), nrow=1)
        L[1,i]<- 1 # select one!!! component of weights vector for a test
        #s.est<- as.numeric((1-cc)*t(mu_est)%*% Q.est %*% mu_est-cc)
        eta.est<-((s.est+cc)/s.est)*L%*%Q.est%*%mu_est/as.numeric(t(mu_est)%*%Q.est%*%mu_est)
        w.est<-(L%*%iS%*%ones)/sum(tones%*%iS%*%ones)+gamma^{-1}*s.est*eta.est # weight (component of weights vector)

        ############# BDOPS, formula 15

        Omega.Lest<- (((1-cc)/(s.est+cc)+ (s.est+cc)/gamma)/gamma + V.est)*(1-cc)*L%*%Q.est%*%t(L)+
          gamma^{-2}*(2*(1-cc)*cc^3/(s.est+cc)^2+ 4*(1-cc)*cc*s.est*(s.est+2*cc)/(s.est+cc)^2 +
                      2*(1-cc)*cc^2*(s.est+cc)^2/(s.est^2)-s.est^2)*eta.est%*%t(eta.est)

        TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est)
        low_bound <- w.est - qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
        upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)

        p_value<- pchisq(TDn, df=1, lower.tail=FALSE)

        #return(c(w.est, low_bound, upp_bound, TDn,  p_value)) #if you want to have the components of the weights vector in F. ()
        return(c(w_EU_shr[i], low_bound, upp_bound, TDn,  p_value)) #first column= shrinkage estimator for EU Portfolio weights
    }

  return(T.dens)
}

stopCluster(cl)


######################################################################################################

 # Example of input
 n<-3e2 # number of realizations
 p<-.5*n # number of assets
 b<-rep(1/p,p)
 gamma<-1
 x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)


#### The function ##############

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
 V_hat_c <- hdsp:::V_hat_c_fast(ones=ones, invSS=invSS, tones=tones, c=cc)  # V_hat_c_fast is V.est
 Q_n_hat <- hdsp:::Q_hat_n_fast(invSS=invSS, Ip=ones, tIp=tones)  # Q_n_hat is Q.est
 s_hat_c <- as.numeric((1-cc)*t(mu_est)%*% Q_n_hat %*% mu_est-cc) # s.est
 R_hat_GMV<-(tones%*% invSS %*% mu_est)/sum(tones%*%invSS%*%ones) # R.est, returns of EU portfolio

 ######## calculate shrinkage weights for EU or GMVP ##############
 V_hat_b <- hdsp:::V_b(Sigma=cov_mtrx, b=b)  # Vb.est
 R_hat_b <- hdsp:::R_b(mu=mu_est, b=b) # Rb.est


 W_EU_hat <- as.vector(
   (invSS %*% ones)/as.numeric(tones %*% invSS %*% ones) +
     Q_n_hat %*% mu_est/gamma,
   mode = 'numeric')

 # alpha_EU
 al <- hdsp:::alpha_hat_star_c_fast(gamma=gamma, c=cc, s=s_hat_c, b=b, R_GMV=R_hat_GMV,
                                    R_b=R_hat_b, V_c=V_hat_c, V_b=V_hat_b)
 weights <- al*W_EU_hat + (1-al)*b # w_EU_shr



 #### Loop ####

 T.dens = foreach(i=1:p, .combine="rbind", .inorder=FALSE, .packages=c("MASS")) %dopar%
   {
     ## EU case
     L <- matrix(rep(0,p), nrow=1)
     L[1,i]<- 1 # select one!!! component of weights vector for a test
     #s.est<- as.numeric((1-cc)*t(mu_est)%*% Q.est %*% mu_est-cc)
     eta.est<-((s.est+cc)/s.est)*L%*%Q.est%*%mu_est/as.numeric(t(mu_est)%*%Q.est%*%mu_est)
     w.est<-(L%*%iS%*%ones)/sum(tones%*%iS%*%ones)+gamma^{-1}*s.est*eta.est # weight (component of weights vector)

     ############# BDOPS, formula 15

     Omega.Lest<- (((1-cc)/(s.est+cc)+ (s.est+cc)/gamma)/gamma + V.est)*(1-cc)*L%*%Q.est%*%t(L)+
       gamma^{-2}*(2*(1-cc)*cc^3/(s.est+cc)^2+ 4*(1-cc)*cc*s.est*(s.est+2*cc)/(s.est+cc)^2 +
                     2*(1-cc)*cc^2*(s.est+cc)^2/(s.est^2)-s.est^2)*eta.est%*%t(eta.est)

     TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est)
     low_bound <- w.est - qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
     upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)

     p_value<- pchisq(TDn, df=1, lower.tail=FALSE)

     #return(c(w.est, low_bound, upp_bound, TDn,  p_value)) #if you want to have the components of the weights vector in F. ()
     return(c(w_EU_shr[i], low_bound, upp_bound, TDn,  p_value)) #first column= shrinkage estimator for EU Portfolio weights
   }






