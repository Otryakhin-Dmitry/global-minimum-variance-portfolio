#inference for each asset (BDPT Sampling distributions ...2020)

test_each_asset <- function(gamma, X_asset, alph){

  library("doParallel")
  no_cores <- detectCores()
  library("foreach")
  cl<-makeCluster(no_cores)

  registerDoParallel(cl)

  p <- ncol(X_asset)
  n <- nrow(X_asset)
  cc<- p/n
  b <-rep(1/p, p) #for shrinkage, target

  ones<-matrix(1,p,1)
  tones<-t(ones)

  T.dens = foreach(i=1:p, .combine="rbind", .inorder=FALSE, .packages=c("MASS")) %dopar%
    {
      L <- matrix(rep(0,p), nrow=1)
      L[1,i]<- 1 # select one!!! component of weights vector for a test
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
      ####  for calculating shrinkage GMVP weights (can be remove if you import your own function)
      w_GMVP_whole<- iS%*%ones/sum(tones%*%iS%*%ones) # sample estimator for GMVP weights
      Lb<- (1-cc)*t(b)%*%S%*%b*sum(tones%*%iS%*%ones) -1 #relative loss of GMVP f.(16) IEEE 2019
      alpha_GMVP<- as.numeric((1-cc)*Lb/(cc+(1-cc)*Lb)) #shrinkage intensity f.(16) IEEE 2019
      w_GMV_shr <- alpha_GMVP*w_GMVP_whole + (1-alpha_GMVP)*b # f.(17) IEEE 2019
      #end
      #######################

      #########################
      ####  for calculating shrinkage EU weights (can be remove if you import your own function)
      w_EU_whole<- iS%*%ones/sum(tones%*%iS%*%ones)+gamma^{-1}*Q.est%*%mu_est # sample estimator for EU weights f.(6) IEEE 2021
      alpha_EU<- as.numeric(gamma^{-1}*((R.est-Rb.est)*(1+1/(1-cc))+ gamma*(Vb.est-V.est)+gamma^{-1}*s.est/(1-cc))/(V.est/(1-cc)-2*(V.est+gamma^{-1}*(Rb.est-R.est)/(1-cc))+ gamma^{-2}*(s.est/(1-cc)^3+cc/(1-cc)^3)+Vb.est))
      #f.(25) IEEE 21
      w_EU_shr <- alpha_EU*w_EU_whole + (1-alpha_EU)*b # f.(26) IEEE 2021
      #end
      #######################


      if (gamma == Inf){ # for GMVP

        w.est <- (L%*%iS%*%ones)/sum(tones%*%iS%*%ones) #(16) for GMVP
        ############# BDOPS, formula 15
        Omega.Lest <- V.est*(1-cc)*L%*%Q.est%*%t(L) # (17) IEEE 2021 simplifies to that)

        TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est) #test statistics

        low_bound <- w.est -qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
        upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)

        p_value <- pchisq(TDn, df=1, lower.tail=FALSE)

        #return(c(w.est, low_bound, upp_bound, TDn,  p_value)) #if you want to have the components of the weights vector in F. ()
        return(c(w_GMV_shr[i], low_bound, upp_bound, TDn,  p_value))

        }

      else { ## EU case
        #s.est<- as.numeric((1-cc)*t(mu_est)%*% Q.est %*% mu_est-cc)
        eta.est<-((s.est+cc)/s.est)*L%*%Q.est%*%mu_est/as.numeric(t(mu_est)%*%Q.est%*%mu_est)
        w.est<-(L%*%iS%*%ones)/sum(tones%*%iS%*%ones)+gamma^{-1}*s.est*eta.est # weight (component of weights vector)

       ############# BDOPS, formula 15

        Omega.Lest<- (((1-cc)/(s.est+cc)+ (s.est+cc)*gamma^{-1})*gamma^{-1}+
                        V.est)*(1-cc)*L%*%Q.est%*%t(L)+
          gamma^{-2}*(2*(1-cc)*cc^3/(s.est+cc)^2+ 4*(1-cc)*cc*s.est*(s.est+2*cc)/(s.est+cc)^2+2*(1-cc)*cc^2*(s.est+cc)^2/(s.est^2)-s.est^2)*eta.est%*%t(eta.est)

        TDn<- (n-p)*t(w.est)%*%solve(Omega.Lest)%*%(w.est)
        low_bound <- w.est - qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)
        upp_bound <- w.est + qnorm(1-alph/2)*sqrt(Omega.Lest)/sqrt(n-p)

        p_value<- pchisq(TDn, df=1, lower.tail=FALSE)

        #return(c(w.est, low_bound, upp_bound, TDn,  p_value)) #if you want to have the components of the weights vector in F. ()
        return(c(w_EU_shr[i], low_bound, upp_bound, TDn,  p_value)) #first column= shrinkage estimator for EU Portfolio weights
      }
    }

  stopCluster(cl)
  return(T.dens)
}

###### not for the package, only for check ######
###### further--> check the performance of the above function, simulate entry matrix
library(MASS)

p <- 300
cc<-0.8
n<- round(p/cc)

mu <- seq(0.2,-0.2, length.out=p) #simulate mu
#####---Covariance matrix from Wishart distr. with given eigenvalues---#####
q<-7.65 # then conditionindex = approx. 450
lambda<- rep(NA, p) #simulate some eigenvalues
for(i in 1:p){
  lambda[i]<- 0.1*exp(q*cc*(i-1)/p)
}
Eigen.matrix <- diag(lambda)

sq.Eigen.matrix <-sqrt(Eigen.matrix)
###---Wishart distribution----####
Z0<-matrix(rnorm(p*n*n),p,n*n)
Wish.matr <- Z0%*%t(Z0)
####----spectral decomposition----####
U <- eigen(Wish.matr)$vectors
#####---Covariance matrix from Wishart distr. with given eigenvalues---#####
cov.mat.H0 <- U%*% Eigen.matrix %*%t(U)

#set.seed(3)
X <- mvrnorm(n,mu, cov.mat.H0)#<--- simulated matrix of asset returns

a<-test_each_asset(5, X, 0.05)#<--- call the function
a[1:15, ]
#sum(a[,1]) # weights should sum up to 1
#which(a[,1] > 0.05)

