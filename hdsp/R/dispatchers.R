# EU portfolio dispatcher

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





