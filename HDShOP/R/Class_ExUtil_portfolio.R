#' S3 class MeanVar_portfolio
#'
#'
#' Class MeanVar_portfolio is designed to construct mean-variance portfolios
#' with provided estimators of the mean vector, covariance matrix, and inverse covariance matrix.
#' It includes the following elements:
#'
#' @section Slots:
#'
#' | Element | Description |
#' | --- | --- |
#' | call | the function call with which it was created |
#' | cov_mtrx | the sample covariance matrix of the asset returns |
#' | inv_cov_mtrx | the inverse of the sample covariance matrix |
#' | means | sample mean vector of the asset returns |
#' | weights | portfolio weights |
#' | Port_Var | portfolio variance |
#' | Port_mean_return | expected portfolio return |
#' | Sharpe | portfolio Sharpe ratio |
#' @md
#'
#' @seealso [summary.MeanVar_portfolio] summary method for the class,
#' [new_MeanVar_portfolio] class constructor, [validate_MeanVar_portfolio] class validator,
#' [MeanVar_portfolio] class helper.
#' @name Class_MeanVar_portfolio
NULL
