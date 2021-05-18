#' A set of tools for shrinkage estimators and Expected Utility portfolios
#'
#' Package HDShOP has the following three important functions: \code{\link{EUShrinkPortfolio}},
#' \code{\link{CovarEstim}} and \code{\link{MeanEstim}}. EUShrinkPortfolio creates Expected Utility
#' portfolios using shrinkage estimation methods including ones for mean vectors, direct and inverse
#' covariance matrices and portfolio weights. CovarEstim computes sample covariance matrices while MeanEstim-
#' mean vectors. Every of these three functions is supplied a name of the method to use to perform estimation.
#' All portfolios are stored in objects of class ExUtil_portfolio. Additionally, most portfolios have
#' a subclass, specific to their kind, that inherits from ExUtil_portfolio. For the later class, constructor,
#' validator and helper functions are available, so that custom EU portfolios may be coded by users.
#'
#'
#' @section Methods:
#'
#' MeanEstim: \insertCite{BOP2019}{HDShOP}, James-Stein and Bayes-Stein estimators \insertCite{Jorion1986}{HDShOP}.
#'
#' CovarEstim: \insertCite{BGP2014}{HDShOP}, \insertCite{LW2020}{HDShOP}.
#'
#' EUShrinkPortfolio: all above plus \insertCite{BDOPS2021}{HDShOP}, \insertCite{BDPS2019}{HDShOP}
#' and \insertCite{BGP2016}{HDShOP}.
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
