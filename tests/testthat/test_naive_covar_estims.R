# context("The simplest covariance estimator")

###############################
p <- 5 # number of assets
n <- 1e6 # number of realizations
###############################


#### Test with zero mean

x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)
SS <- Sigma_sample_estimator(x)


test_that("Non-diagonal elements are small when components are independent", {
  expect_lt(max(abs(SS[2,3]), abs(SS[2,4]), abs(SS[2,5])), 0.05)
})

test_that("Diagonal elements are close to 1 when components are independent", {
  expect_lt(1-min(abs(SS[1,1]), abs(SS[2,2]), abs(SS[3,3]),
                  abs(SS[4,4]), abs(SS[5,5])), 0.05)
})



#### Test with non-zero mean

x <-matrix(data = rnorm(n*p, mean=10), nrow = p, ncol = n)
SS <- Sigma_sample_estimator(x)

test_that("Non-diagonal elements are small when components are independent", {
  expect_lt(max(abs(SS[2,3]), abs(SS[2,4]), abs(SS[2,5])), 0.05)
})

test_that("Diagonal elements are close to 1 when components are independent", {
  expect_lt(1-min(abs(SS[1,1]), abs(SS[2,2]), abs(SS[3,3]),
                  abs(SS[4,4]), abs(SS[5,5])), 0.05)
})


# Test with cov.shrink

test_that("cov.shrink produces about the same matrix as Sigma_sample_estimator", {

  if (!requireNamespace("corpcor", quietly =TRUE)) skip("corpcor is not installed")
  library(corpcor)

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  ####
  n <- 1e3 # number of realizations
  p <- 0.8*n # number of assets

  x <- matrix(data = rnorm(n*p, mean=10), nrow = p, ncol = n)

  s2 <- cov.shrink(t(x), lambda=0, verbose=FALSE)
  s1 <- Sigma_sample_estimator(x)
  Expr <- abs(s2[1:6,1:6] - s1[1:6,1:6]) / abs(s2[1:6,1:6])<0.2

  expect_equal(sum(Expr), 36)
})




