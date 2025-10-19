# context("BDPS19 and BDOPS21 return equivalent portfolios when gamma==Inf")

library(MASS)
n<-3e2 # number of realizations
p<-.5*n # number of assets
b<-rep(1/p,p)


Mtrx <- RandCovMtrx(p=p)
x <- t(mvrnorm(n=n , mu=rep(0,p), Sigma=Mtrx))

test_GMV <- new_GMV_portfolio_weights_BDPS19(x=x, b=b, beta=0.05)
test_MV <- new_MV_portfolio_weights_BDOPS21(x=x, gamma=Inf, b=b, beta=0.05)

test_that("Elements of outputs GMV and EU portfolios are equal", {

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  # Weights
  expect_equal(test_GMV$weights, test_MV$weights, ignore_attr = TRUE)
  expect_equal(test_GMV$means, test_MV$means, ignore_attr = TRUE)
  expect_equal(test_GMV$cov_mtrx, test_MV$cov_mtrx, ignore_attr = TRUE)
  expect_equal(test_GMV$inv_cov_mtrx, test_MV$inv_cov_mtrx, ignore_attr = TRUE)
})
