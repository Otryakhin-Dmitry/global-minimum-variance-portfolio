



#### EU portfolio dispatcher

#' Shrinkage portfolio
#'
#' Function dispatcher for portfolio construction
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='weights', b=b)
#' str(test)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='Bayes-Stein')
#' str(test)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='James-Stein')
#' str(test)
EUShrinkPortfolio <- function(x, gamma, type, subtype, ...) {

  if(!is.numeric(gamma) || is.na(gamma)) stop("gamma is not numeric")
  if(gamma==Inf) stop("GMVP methods are absent")

  if(type=='mean') {
    if(subtype=='James-Stein') {
      output <- new_ExUtil_portfolio_mean_JamesStein(x=x, gamma=gamma, ...)
    }
    if(subtype=='Bayes-Stein') {
      output <- new_ExUtil_portfolio_mean_BayesStein(x=x, gamma=gamma, ...)
    }
  }

  if(type=='cov') {
    if(subtype=='LW02') {
      output <- new_ExUtil_portfolio_cov_LW02(x=x, gamma=gamma, ...)
    }
    if(subtype=='BGP14') {
      output <- new_ExUtil_portfolio_cov_BGP14(x=x, gamma=gamma, ...)
    }
  }

  if(type=='inv_cov') {
    if(subtype=='BGP16') {
      output <- new_ExUtil_portfolio_pm_BGP16(x=x, gamma=gamma, ...)
    }
  }

  if(type=='weights') {
    output <- new_ExUtil_portfolio(x=x, gamma=gamma, ...)
  }
  output
}



#### Covariance shrinkage ####

#' Covariance matrix estimator
#'
#' Function dispatcher for portfolio construction
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' Mtrx_naive <- CovarEstim(x, subtype="naive")
#'
#' TM <- matrix(0, p, p)
#' diag(TM) <- 1
#' Mtrx_bgp <- CovarEstim(x, subtype="BGP14", TM=TM)
#'
#' Mtrx_bgp <- CovarEstim(x, subtype="LW20", TM=TM)
CovarEstim <- function(x, subtype, ...)
{
    if(subtype=='naive') {
      output <- Sigma_sample_estimator(x=x)
    }
    if(subtype=='BGP14') {
      SCM <- Sigma_sample_estimator(x=x)
      output <- CovShrinkBGP14(n=n, TM=TM, SCM=SCM)
    }
    if(subtype=='LW20') {
      output <- nonlin_shrinkLW(x=x)
    }
    output
}



#### Mean vector shrinkage ####

#' Mean vector shrinkage estimator
#'
#' Function dispatcher for portfolio construction
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#'
#' Mean_naive <- MeanEstim(x, subtype="naive")
#' Mean_BGP <- MeanEstim(x, subtype="naive")
#'
#' mu_0 <- rep(1/p, p)
#'
MeanEstim <- function(x, subtype, ...)
{
  if(subtype=='naive') {
    n <- ncol(x)
    p <- nrow(x)
    output <- .rowMeans(x=x, m=p, n=n)
  }
  if(subtype=='BOP19') {
    output <- mean_bop19(x, mu_0)
  }
  if(subtype=='js') {
    output <- mean_js(x=x, mu_0=mu_0)
  }
  if(subtype=='bs') {
    output <- mean_bs(x=x, mu_0=mu_0)
  }
  output
}





