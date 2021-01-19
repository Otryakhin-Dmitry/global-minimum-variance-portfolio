#' A set of tools for shrinkage estimators and Expected Utility portfolios
#'
#' Package hdsp has the following three important functions: \code{\link{EUShrinkPortfolio}},
#' \code{\link{CovarEstim}} and \code{\link{MeanEstim}}. EUShrinkPortfolio creates Expected Utility
#' portfolios using shrinkage estimation methods including ones for mean vectors, direct and inverse
#' covariance matrices and portfolio weights. CovarEstim computes sample covariance matrices while MeanEstim-
#' mean vectors. Every of these three functions is supplied a name of the method to use to perform estimation.
#' Portfolios are stored in objects of class ExUtil_portfolio. For this class constructor, validator and
#' helper functions are available.
#'
#'
#' @section Methods:
#'
#' MeanEstim: \insertCite{BOP2019}{hdsp}, James-Stein and Bayes-Stein estimators \insertCite{Jorion1986}{hdsp}.
#'
#' CovarEstim: \insertCite{BGP2014}{hdsp}, \insertCite{LW2020}{hdsp}.
#'
#' EUShrinkPortfolio: all above plus \insertCite{BDOPS2020}{hdsp} and \insertCite{BGP2016}{hdsp}.
#'
#'
#' @name hdsp-package
#' @docType package
#' @references
#' \insertAllCited{}
NULL

# Importing from Rdpack (for working with references)
# to avoid getting a warning from "R CMD check"
#' @importFrom Rdpack reprompt
NULL
