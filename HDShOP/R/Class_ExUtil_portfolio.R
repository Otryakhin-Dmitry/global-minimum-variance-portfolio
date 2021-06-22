#' S3 class ExUtil_portfolio
#'
#'
#' Class ExUtil_portfolio is designed to describe Expected Utility portfolios with
#' existing direct and inverse covariance matrices. It includes the following elements:
#'
#' @section Slots:
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the asset returns |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean vector estimate for the asset returns |
#' | weights | portfolio weights |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | expected portfolio return |
#' | Sharpe | portfolio Sharpe ratio |
#' @md
#'
#' @seealso [summary.ExUtil_portfolio] summary method for the class,
#' [new_MV_portfolio_custom] class constructor, [validate_MV_portfolio] class validator,
#' [MV_portfolio_custom] class helper.
#' @name ExUtil_portfolio
NULL
