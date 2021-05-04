

#### EU portfolio dispatcher

#' Shrinkage Expected Utility portfolio
#'
#' The main function for EU portfolio construction. It is a dispatcher using methods according
#' to arguments type and subtype.
#'
#' The available estimation methods are:
#'
#' | Function | Paper | Type | Subtype |
#' | --- | --- | --- | --- |
#' | \code{\link{new_ExUtil_portfolio_mean_BayesStein}} | Jorion 1986 | mean | Bayes-Stein |
#' | \code{\link{new_ExUtil_portfolio_mean_JamesStein}} | Jorion 1986 | mean | James-Stein |
#' | \code{\link{new_ExUtil_portfolio_mean_BOP19}} | Bodnar et al 2019 | mean | BOP19 |
#' | \code{\link{new_ExUtil_portfolio_cov_LW02}} | Ledoit & Wolf 2020 | cov | LW02 |
#' | \code{\link{new_ExUtil_portfolio_cov_BGP14}} | Bodnar et al 2014 | cov | BGP14 |
#' | \code{\link{new_ExUtil_portfolio_icov_BGP16}} | Bodnar et al 2016 | inv_cov | BGP16 |
#' | \code{\link{new_ExUtil_portfolio_weights_BDOPS20}} | Bodnar et al 2020 | weights | |
#'
#' @md
#' @param x a matrix or a data frame of asset returns. Rows represent different assets, columns- observations.
#' @param gamma a numeric variable. Investors attitude towards risk.
#' @param type a character. The type of methods to use to construct the portfolio.
#' @param subtype a character. The exact method to use within the type.
#' @param ... arguments to pass to portfolio constructors
#'
#' @return A portfolio in the form of an object of class ExUtil_portfolio potentially with a subclass.
#' See \code{\link{new_ExUtil_portfolio_custom}} for the details of the class.
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='weights', b=b, alph = 0.05)
#' str(test)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='Bayes-Stein')
#' str(test)
#'
#' test <- EUShrinkPortfolio(x=x, gamma=gamma, type='mean', subtype='James-Stein')
#' str(test)
#' @export
EUShrinkPortfolio <- function(x, gamma, type, subtype, ...) {

  if(!is.numeric(gamma) || is.na(gamma)) stop("gamma is not numeric")
  if(gamma==Inf) stop("GMVP methods are absent")

  if(!is.character(type)) stop("type is not character")

  if(type=='mean') {

    if(subtype=='BOP19') {
      output <- new_ExUtil_portfolio_mean_BOP19(x=x, gamma=gamma, ...)
      return(output)

    } else if(subtype=='James-Stein') {
      output <- new_ExUtil_portfolio_mean_JamesStein(x=x, gamma=gamma, ...)
      return(output)

    } else if(subtype=='Bayes-Stein') {
      output <- new_ExUtil_portfolio_mean_BayesStein(x=x, gamma=gamma, ...)
      return(output)

    } else {
      stop(paste('Invalid subtype for type',type,sep=" "))
    }
  } else if(type=='cov') {

    if(subtype=='LW02') {
      output <- new_ExUtil_portfolio_cov_LW02(x=x, gamma=gamma, ...)
      return(output)

    } else if(subtype=='BGP14') {
      output <- new_ExUtil_portfolio_cov_BGP14(x=x, gamma=gamma, ...)
      return(output)

    } else {
      stop(paste('Invalid subtype for type',type,sep=" "))
    }
  } else if(type=='inv_cov') {

    if(subtype=='BGP16') {
      output <- new_ExUtil_portfolio_icov_BGP16(x=x, gamma=gamma, ...)
      return(output)

    } else {
      stop(paste('Invalid subtype for type',type,sep=" "))
    }
  } else  if(type=='weights') {

    output <- new_ExUtil_portfolio_weights_BDOPS20(x=x, gamma=gamma, ...)
    return(output)
  } else {

    stop(paste('Invalid type:',type,sep=" "))
  }
}



#### Covariance shrinkage ####

#' Covariance matrix estimator
#'
#' Essentially it is a function dispatcher for covariance estimation that chooses a method according to the type argument.
#'
#' The available estimation methods are:
#'
#' | Function | Paper | Type |
#' | --- | --- | --- |
#' | \code{\link{Sigma_sample_estimator}} |  | naive |
#' | \code{\link{CovShrinkBGP14}} | Bodnar et al 2014 | BGP14 |
#' | \code{\link{nonlin_shrinkLW}} | Ledoit & Wolf 2020| LW20 |
#' @md
#'
#' @param x a matrix of asset returns. Rows represent different assets, columns- observations.
#' @param type a character. The the estimation method to use.
#' @param ... arguments to pass to estimators
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' Mtrx_naive <- CovarEstim(x, type="naive")
#'
#' TM <- matrix(0, p, p)
#' diag(TM) <- 1
#' Mtrx_bgp <- CovarEstim(x, type="BGP14", TM=TM, SCM=Mtrx_naive)
#'
#' Mtrx_bgp <- CovarEstim(x, type="LW20")
#' @export
CovarEstim <- function(x, type, ...)
{
    if(type=='naive') {
      output <- Sigma_sample_estimator(x=x)
    }
    if(type=='BGP14') {
      SCM <- Sigma_sample_estimator(x=x)
      n<-ncol(x)
      output <- CovShrinkBGP14(n, ...)
    }
    if(type=='LW20') {
      output <- nonlin_shrinkLW(x=x)
    }
    output
}



#### Mean vector shrinkage ####

#' Mean vector shrinkage estimator
#'
#' Essentially it is a function dispatcher for mean value estimators that chooses a method according to the type argument.
#'
#' The available estimation methods are:
#'
#' | Function | Paper | Type |
#' | --- | --- | --- |
#' | .rowMeans |  | naive |
#' | \code{\link{mean_bs}} | Jorion 1986 | bs |
#' | \code{\link{mean_js}} | Jorion 1986 | js |
#' | \code{\link{mean_bop19}} | Bodnar et al 2019 | BOP19 |
#' @md
#'
#' @inheritParams CovarEstim
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' Mean_naive <- MeanEstim(x, type="naive")
#'
#' mu_0 <- rep(1/p, p)
#' Mean_BOP <- MeanEstim(x, type="BOP19", mu_0=mu_0)
#'
#' @export
MeanEstim <- function(x, type, ...)
{
  if(type=='naive') {
    n <- ncol(x)
    p <- nrow(x)
    output <- .rowMeans(x=x, m=p, n=n)
  }
  if(type=='BOP19') {
    output <- mean_bop19(x=x, ...)
  }
  if(type=='js') {
    output <- mean_js(x=x, ...)
  }
  if(type=='bs') {
    output <- mean_bs(x=x, ...)
  }
  output
}





