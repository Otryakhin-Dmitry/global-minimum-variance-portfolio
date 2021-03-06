% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3_tradit_portfol.R
\name{new_MV_portfolio_traditional}
\alias{new_MV_portfolio_traditional}
\title{Traditional mean-variance portfolio}
\usage{
new_MV_portfolio_traditional(x, gamma)
}
\arguments{
\item{x}{a matrix or a data frame of asset returns. Rows represent different assets, columns- observations.}

\item{gamma}{a numeric variable. Coefficient of risk aversion.}
}
\value{
an object of class MeanVar_portfolio\tabular{ll}{
   Element \tab Description \cr
   call \tab the function call with which it was created \cr
   cov_mtrx \tab the sample covariance matrix of asset returns \cr
   inv_cov_mtrx \tab the inverse of the sample covariance matrix \cr
   means \tab sample mean estimate for the asset returns \cr
   W_mv_hat \tab sample estimate of portfolio weights \cr
   Port_Var \tab portfolio variance \cr
   Port_mean_return \tab expected portfolio return \cr
   Sharpe \tab portfolio Sharpe ratio \cr
}
}
\description{
Mean-variance portfolios with the traditional (sample) estimators for the mean
vector and the covariance matrix of asset returns. For more details of the method,
see \code{\link{MVShrinkPortfolio}}.
}
\examples{
n<-3e2 # number of realizations
p<-.5*n # number of assets
gamma<-1

x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)

test <- new_MV_portfolio_traditional(x=x, gamma=gamma)
str(test)
}
