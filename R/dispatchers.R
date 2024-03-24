

#### EU portfolio dispatcher

#' Shrinkage mean-variance portfolio
#'
#' The main function for mean-variance (also known as expected utility)
#' portfolio construction. It is a dispatcher using methods according to
#' argument type, values of gamma and dimensionality of matrix x.
#'
#' The sample estimator of the mean-variance portfolio weights, which results in
#' a traditional mean-variance portfolio, is calculated by
#' \deqn{\hat w_{MV} = \frac{S^{-1} 1}{1' S^{-1} 1} +
#'       \gamma^{-1} \hat Q \bar x,}
#' where \eqn{S^{-1}} and \eqn{\bar x} are the inverse of the sample covariance
#' matrix and the sample mean vector of asset returns respectively, \eqn{\gamma}
#' is the coefficient of risk aversion and \eqn{\hat Q} is given by
#' \deqn{\hat Q = S^{-1} - \frac{S^{-1} 1 1' S^{-1}}{1' S^{-1} 1} .}
#' In the case when \eqn{p>n}, \eqn{S^{-1}} becomes \eqn{S^{+}}- Moore-Penrose
#' inverse. The shrinkage estimator for the mean-variance portfolio weights in
#' a high-dimensional setting is given by
#' \deqn{\hat w_{ShMV} = \hat \alpha \hat w_{MV} + (1- \hat \alpha)b,}
#' where \eqn{\hat \alpha} is the estimated shrinkage intensity and \eqn{b} is
#' a target vector with the sum of the elements equal to one.
#'
#' In the case \eqn{\gamma \neq \infty}, \eqn{\hat{\alpha}} is computed following
#' Eq. (2.22) of \insertCite{BOP16;textual}{HDShOP} for c<1 and following
#' Eq. (2.29) of \insertCite{BOP16;textual}{HDShOP} for c>1.
#'
#' The case of a fully risk averse investor (\eqn{\gamma=\infty}) leads to the
#' traditional global minimum variance (GMV) portfolio with the weights given by
#' \deqn{\hat w_{GMV} = \frac{S^{-1} 1}{1' S^{-1} 1} .}
#' The shrinkage estimator for the GMV portfolio is then calculated by
#' \deqn{\hat w_{ShGMV} = \hat\alpha \hat w_{GMV} + (1-\hat \alpha)b,}
#' with \eqn{\hat{\alpha}}  given in
#' Eq. (2.31) of \insertCite{BPS2018;textual}{HDShOP} for c<1 and in
#' Eq. (2.33) of \insertCite{BPS2018;textual}{HDShOP} for c>1.
#'
#' These estimation methods are available as separate functions employed by
#' MVShrinkPortfolio accordingly to the following parameter configurations:
#'
#' | Function | Paper | Type | gamma | p/n |
#' | --- | --- | --- | --- | --- |
#' | \code{\link{new_MV_portfolio_weights_BDOPS21}} | \insertCite{BOP16;textual}{HDShOP} | shrinkage | < Inf | <1 |
#' | \code{\link{new_MV_portfolio_weights_BDOPS21_pgn}} | \insertCite{BOP16;textual}{HDShOP} | shrinkage | < Inf | >1 |
#' | \code{\link{new_GMV_portfolio_weights_BDPS19}} | \insertCite{BPS2018;textual}{HDShOP} | shrinkage | Inf | <1 |
#' | \code{\link{new_GMV_portfolio_weights_BDPS19_pgn}} | \insertCite{BPS2018;textual}{HDShOP} | shrinkage | Inf | >1 |
#' | \code{\link{new_MV_portfolio_traditional}} |  | traditional | > 0 | <1 |
#' | \code{\link{new_MV_portfolio_traditional_pgn}} |  | traditional | > 0 | >1 |
#' @md
#' @param x a p by n matrix or a data frame of asset returns. Rows represent
#' different assets, columns -- observations.
#' @param gamma a numeric variable. Coefficient of risk aversion.
#' @param type a character. The type of methods to use to construct the
#' portfolio.
#' @param ... arguments to pass to portfolio constructors
#'
#' @return A portfolio in the form of an object of class MeanVar_portfolio
#' potentially with a subclass.
#' See \code{\link{Class_MeanVar_portfolio}} for the details of the class.
#' @references \insertAllCited{}
#' @examples
#' n<-3e2 # number of realizations
#' gamma<-1
#'
#' # The case p<n
#'
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=gamma,
#'                           type='shrinkage', b=b, beta = 0.05)
#' str(test)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=Inf,
#'                           type='shrinkage', b=b, beta = 0.05)
#' str(test)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=gamma, type='traditional')
#' str(test)
#'
#' # The case p>n
#'
#' p<-1.2*n # Re-define the number of assets
#' b<-rep(1/p,p)
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=gamma, type='shrinkage',
#'                           b=b, beta = 0.05)
#' str(test)
#'
#' test <- MVShrinkPortfolio(x=x, gamma=Inf, type='shrinkage',
#'                           b=b, beta = 0.05)
#' str(test)
#'
#' @export
MVShrinkPortfolio <- function(x, gamma,
                              type=c('shrinkage', 'traditional'), ...) {

  if(!is.numeric(gamma) || is.na(gamma)) stop("gamma is not numeric")
  if(!is.character(type)) stop("type is not character")

  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.array(x) || length(dim(x)) != 2L)
    stop("'x' must be a matrix or a dataframe")

  cl <- match.call()

  if(type=='traditional') {

    if(nrow(x) >= ncol(x)) {
      output <- new_MV_portfolio_traditional_pgn(x=x, gamma=gamma)
    } else {
      output <- new_MV_portfolio_traditional(x=x, gamma=gamma)
    }


  } else  if(type=='shrinkage') {

    if(gamma != Inf) {

      if(nrow(x) < ncol(x)) {
        output <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=gamma, ...)
      } else if(nrow(x) > ncol(x)) {
        output <- new_MV_portfolio_weights_BDOPS21_pgn(x=x, gamma=gamma, ...)
      } else if(nrow(x) == ncol(x)) stop('do not know how to treat a square x')

    } else {

      if(nrow(x) < ncol(x)) {
        output <- new_GMV_portfolio_weights_BDPS19(x=x, ...)
      } else if(nrow(x) > ncol(x)) {
        output <- new_GMV_portfolio_weights_BDPS19_pgn(x=x, ...)
      } else if(nrow(x) == ncol(x)) stop('do not know how to treat a square x')

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
#' @param x a p by n matrix or a data frame of asset returns. Rows represent
#' different assets, columns -- observations.
#' @param type a character. The estimation method to be used.
#' @param ... arguments to pass to estimators
#' @return an object of class matrix
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' Mtrx_trad <- CovarEstim(x, type="trad")
#'
#' TM <- matrix(0, p, p)
#' diag(TM) <- 1
#' Mtrx_bgp <- CovarEstim(x, type="BGP14", TM=TM)
#'
#' Mtrx_lw <- CovarEstim(x, type="LW20")
#' @export
CovarEstim <- function(x, type=c('trad', 'BGP14', 'LW20'), ...)
{
    if (is.data.frame(x)) x <- as.matrix(x)
    if (!is.array(x) || length(dim(x)) != 2L)
      stop("'x' must be a matrix or a dataframe")

    if(type=='trad') {
      output <- Sigma_sample_estimator(x=x)
    }
    if(type=='BGP14') {
      ll <- list(...)
      if(is.null(ll$SCM)) {
        SCM <- Sigma_sample_estimator(x=x)
      } else {
        SCM <- as.matrix(ll$SCM)
      }
      n<-ncol(x)
      output <- CovShrinkBGP14(n=n, SCM=SCM, TM=as.matrix(ll$TM))$S
    }
    if(type=='LW20') {
      output <- nonlin_shrinkLW(x=x)
    }
    output
}



#### Mean vector shrinkage ####

#' Mean vector estimator
#'
#' A user-friendly function for estimation of the mean vector.
#' Essentially, it is a function dispatcher for estimation of the mean vector
#' that chooses a method accordingly to the type argument.
#'
#' The available estimation methods for the mean are:
#'
#' | Function | Paper | Type |
#' | --- | --- | --- |
#' | .rowMeans |  | trad |
#' | \code{\link{mean_bs}} | Jorion 1986 | bs |
#' | \code{\link{mean_js}} | Jorion 1986 | js |
#' | \code{\link{mean_bop19}} | Bodnar et al 2019 | BOP19 |
#' @md
#'
#' @inheritParams CovarEstim
#' @return a numeric vector--- a value of the specified estimator of
#' the mean vector.
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
#' Mean_trad <- MeanEstim(x, type="trad")
#'
#' mu_0 <- rep(1/p, p)
#' Mean_BOP <- MeanEstim(x, type="BOP19", mu_0=mu_0)
#'
#' @export
MeanEstim <- function(x, type, ...)
{
  if (is.data.frame(x)) x <- as.matrix(x)
  if (!is.array(x) || length(dim(x)) != 2L)
    stop("'x' must be a matrix or a dataframe")

  if(type=='trad') {
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





