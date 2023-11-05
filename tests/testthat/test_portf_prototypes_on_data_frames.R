context("Portfolio prototypes on data frames")

#### Parameter setting
c <- 0.4
n <- 7e2
p <- c*n
gamma <- 1

mu_n <- rep(0, p)
Sigma <- matrix(0, p, p)
diag(Sigma) <- 1


#### Generate a sample of asset returns
X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)
X <- as.data.frame(X)


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
