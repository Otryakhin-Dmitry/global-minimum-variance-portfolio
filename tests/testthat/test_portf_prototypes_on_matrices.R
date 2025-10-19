# context("Portfolio prototypes")

#### Parameter setting
c <- 0.4
n <- 7e2
p <- c*n
gamma <- 1

mu_n <- rep(0, p)
Sigma <- matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- solve(Sigma)


#### Generate a sample of asset returns
X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)


# weight shrink portfolios
port_BDOPS21 <- new_MV_portfolio_weights_BDOPS21(x=X, gamma=gamma, b=rep(1/p,p), beta = 0.05)
port_BDPS19 <- new_GMV_portfolio_weights_BDPS19(x=X, b=rep(1/p,p), beta=0.05)


#### tests

test_that("All prototypes produce not NULLs", {

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  # Weights
  expect_false(is.null(port_BDOPS21))
  expect_false(is.null(port_BDPS19))
})


test_that("All prototypes are valid MeanVar_portfolio objects", {

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  # Weights
  expect_s3_class(object=validate_MeanVar_portfolio(port_BDOPS21), class="MeanVar_portfolio", exact = FALSE)
  expect_s3_class(object=validate_MeanVar_portfolio(port_BDPS19), class="MeanVar_portfolio", exact = FALSE)
})


