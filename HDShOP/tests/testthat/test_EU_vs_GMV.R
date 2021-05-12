context("BDPS19 and BDOPS21 return equivalent portfolios when gamma==Inf")

library(MASS)
n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)


Mtrx <- RandCovMtrx(n=n, p=p, q=20.55)
x <- t(mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))

test_GMV <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, alph=0.05)
test_EU <- new_ExUtil_portfolio_weights_BDOPS21(x=x, gamma=Inf, b=b, alph=0.05)

test_that("Elements of outputs GMV and EU portfolios are equal", {

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  # Weights
  expect_equivalent(test_GMV$weights, test_EU$weights)
  expect_equivalent(test_GMV$means, test_EU$means)
  expect_equivalent(test_GMV$cov_mtrx, test_EU$cov_mtrx)
  expect_equivalent(test_GMV$inv_cov_mtrx, test_EU$inv_cov_mtrx)
})
