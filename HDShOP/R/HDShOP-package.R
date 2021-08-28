#' A set of tools for shrinkage estimation of mean-variance optimal portfolios
#'
#' Package HDShOP has the following three important functions: \code{\link{MVShrinkPortfolio}},
#' \code{\link{CovarEstim}} and \code{\link{MeanEstim}}. MVShrinkPortfolio creates mean-variance
#' portfolios using shrinkage estimation methods for portfolio weights. CovarEstim computes several estimators
#' of the covariance matrix, while MeanEstim computes several estimators of the mean vector. Each of these three functions is
#' supplied a name of the method used to perform the estimation. All portfolios are stored in objects of
#' class MeanVar_portfolio and some have a subclass, specific to their kind, that inherits from
#' MeanVar_portfolio. For the latter class constructor, validator and helper functions are available, so
#' that custom mean-variance portfolios may be coded by users.
#'
#'
#' @section Methods:
#'
#' MeanEstim: \insertCite{BOP2019}{HDShOP}, James-Stein and Bayes-Stein estimators \insertCite{Jorion1986}{HDShOP}.
#'
#' CovarEstim: \insertCite{BGP2014}{HDShOP}, \insertCite{LW2020}{HDShOP}.
#'
#' MVShrinkPortfolio: \insertCite{BDOPS2021}{HDShOP}, \insertCite{BDPS2019}{HDShOP}.
#'
#'
#' @name HDShOP-package
#' @docType package
#' @references
#' \insertAllCited{}
NULL

# Importing from Rdpack (for working with references)
# to avoid getting a warning from "R CMD check"
#' @importFrom Rdpack reprompt
NULL

# Imports from stats
#' @importFrom stats cov pchisq qnorm pnorm
NULL
