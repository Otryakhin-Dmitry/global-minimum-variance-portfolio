

#### EU portfolio dispatcher

#' Shrinkage Mean Variance portfolio
#'
#' The main function for Mean Variance (also known as Expected Utility) portfolio construction.
#' It is a dispatcher using methods according to argument type.
#'
#' The sample estimator of the mean-variance portfolio weights, which results in
#' a traditional MV portfolio, is calculated by
#' \deqn{\hat w_{MV} = \frac{S^{-1} 1}{1' S^{-1} 1} + \gamma^{-1} \hat Q \bar x \quad ,}
#' where \eqn{S^{-1}} and \eqn{\bar x} are the inverse of the sample covariance
#' matrix and the sample mean vector of asset returns respectively, \eqn{\gamma}
#' is the coefficient of risk aversion and \eqn{\hat Q} is given by
#' \deqn{\hat Q = S^{-1} - \frac{S^{-1} 1 1' S^{-1}}{1' S^{-1} 1} .}
#' Then the shrinkage estimator for MV portfolio weights in a high-dimensional
#' setting is given by \deqn{\hat w_{GSMV} = \hat \alpha \hat w_{MV} + (1- \hat \alpha)b \quad,}
#' where \eqn{\hat \alpha} is the estimated shrinkage intensity and \eqn{b} is
#' a target vector.
#'
#' In the case \eqn{\gamma \neq \infty}, estimation of \eqn{\alpha} is based on
#' a procedure proposed by Eq.2.28 of \insertCite{BOP16;textual}{HDShOP}. The case of a fully risk averse investor
#' (\eqn{\gamma=\infty}) leads  to the construction of the traditional global minimum
#' variance portfolio (GMVP):
#' \deqn{\hat w_{GMV} = \frac{S^{-1} 1}{1' S^{-1} 1} .}
#' The shrinkage estimator for the GMV portfolio is then calculated by
#' \deqn{\hat w_{GS GMV} = \hat \alpha \hat w_{GMV} + (1-\hat \alpha)b \quad,}
#' with \eqn{\hat \alpha} as in Eq.2.31 \insertCite{BPS2018;textual}{HDShOP}.
#'
#' These three estimation methods are available as separate functions dispatched
#' accordingly to the following parameter configurations:
#'
#' | Function | Paper | Type | gamma |
#' | --- | --- | --- | --- |
#' | \code{\link{new_MV_portfolio_weights_BDOPS21}} | Bodnar et al 2021 | shrinkage | < Inf |
#' | \code{\link{new_GMV_portfolio_weights_BDPS19}} | Bodnar et al 2019 | shrinkage | Inf |
#' | \code{\link{new_MV_portfolio_traditional}} |  | traditional | > 0 |
#' @md
#' @param x a matrix or a data frame of asset returns. Rows represent different assets, columns- observations.
#' @param gamma a numeric variable. Coefficient of risk aversion.
#' @param type a character. The type of methods to use to construct the portfolio.
#' @param ... arguments to pass to portfolio constructors
#'
#' @return A portfolio in the form of an object of class MeanVar_portfolio potentially with a subclass.
#' See \code{\link{new_MeanVar_portfolio}} for the details of the class.
#' @references \insertAllCited{}
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=gamma, type='shrinkage', b=b, beta = 0.05)
#' str(test)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=Inf, type='shrinkage', b=b, beta = 0.05)
#' str(test)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=gamma, type='traditional')
#' str(test)
#'
#' @export
MVShrinkPortfolio <- function(x, gamma, type='shrinkage', ...) {

  if(!is.numeric(gamma) || is.na(gamma)) stop("gamma is not numeric")
  if(!is.character(type)) stop("type is not character")

  cl <- match.call()

  if(type=='traditional') {
    output <- new_MV_portfolio_traditional(x=x, gamma=gamma)

  } else  if(type=='shrinkage') {

    if(gamma != Inf) {
      output <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=gamma, ...)

    } else {
      output <- new_GMV_portfolio_weights_BDPS19(x=x, ...)
    }

  } else {

    stop(paste('Invalid type:',type,sep=" "))
  }

  output$call <- cl
  return(output)
}


#### Covariance shrinkage ####

#' Covariance matrix estimator
#'
#' It is a function dispatcher for covariance matrix estimation. One can choose
#' between traditional and shrinkage-based estimators.
#'
#' The available estimation methods are:
#'
#' | Function | Paper | Type |
#' | --- | --- | --- |
#' | \code{\link{Sigma_sample_estimator}} |  | traditional |
#' | \code{\link{CovShrinkBGP14}} | Bodnar et al 2014 | BGP14 |
#' | \code{\link{nonlin_shrinkLW}} | Ledoit & Wolf 2020| LW20 |
#'
#' @md
#'
#' @param x a p by n data matrix. Rows represent different variables, columns- observations.
#' @param type a character. The estimation method to be used.
#' @param ... arguments to pass to estimators
#' @return an object of class matrix
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
#' Mtrx_bgp <- CovarEstim(x, type="BGP14", TM=TM)
#'
#' Mtrx_lw <- CovarEstim(x, type="LW20")
#' @export
CovarEstim <- function(x, type, ...)
{
    if(type=='naive') {
      output <- Sigma_sample_estimator(x=x)
    }
    if(type=='BGP14') {
      ll <- list(...)
      SCM <- Sigma_sample_estimator(x=x)
      n<-ncol(x)
      output <- CovShrinkBGP14(n=n, SCM=SCM, TM=ll$TM)$S
    }
    if(type=='LW20') {
      output <- nonlin_shrinkLW(x=x)
    }
    output
}



#### Mean vector shrinkage ####

#' Mean vector shrinkage estimator
#'
#' A user-friendly function for estimation of the mean vector from a set of vectors.
#' Essentially, it is a function dispatcher for estimation of the mean vector that
#' chooses a method accordingly to the type argument.
#'
#' The available estimation methods for the mean are:
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
#' @return a numeric vector containing the specified estimation of the mean vector.
#' @references
#' \insertRef{Jorion1986}{HDShOP}
#'
#' \insertRef{BOP2019}{HDShOP}
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
    output <- mean_bop19(x=x, ...)$means
  }
  if(type=='js') {
    output <- mean_js(x=x, ...)$means
  }
  if(type=='bs') {
    output <- mean_bs(x=x, ...)$means
  }
  output
}





