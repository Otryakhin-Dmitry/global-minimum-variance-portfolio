context("Portfolio prototypes")

#### Parameter setting
c <- 0.4
n <- 7e2
p <- c*n
gamma <- 1

mu_n<-rep(0, p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- solve(Sigma)

#### Generate a sample of asset returns
X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)

# cov / icov shrink portfolios
port_LW02 <- new_ExUtil_portfolio_cov_LW02(x=X, gamma)
port_BGP14 <-new_ExUtil_portfolio_cov_BGP14(x=X, gamma, TM=Sigma)
port_BGP16 <-new_ExUtil_portfolio_icov_BGP16(x=X, gamma, TM=Sigma)

# mean shrink portfolios
port_BS <- new_ExUtil_portfolio_mean_BayesStein(x=X, gamma=gamma, mu_0=rep(1,p))
port_JS <- new_ExUtil_portfolio_mean_JamesStein(x=X, gamma=gamma, mu_0=rep(1,p))
port_BOP<- new_ExUtil_portfolio_mean_BOP(x=X, gamma=gamma, mu_0=rep(1,p))

# weight shrink portfolios
port_BDOPS20 <- new_ExUtil_portfolio_weights_BDOPS20(x=X, gamma=gamma, b=rep(1/p,p))

#### tests

test_that("All prototypes produce not NULLs", {

  # (i-)Covariance
  expect_false(is.null(port_LW02))
  expect_false(is.null(port_BGP14))
  expect_false(is.null(port_BGP16))

  # Mean
  expect_false(is.null(port_BS))
  expect_false(is.null(port_JS))
  expect_false(is.null(port_BOP))

  # Weights
  expect_false(is.null(port_BDOPS20))
})
