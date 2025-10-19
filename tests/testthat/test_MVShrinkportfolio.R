# context("Portfolio prototypes")


#### Parameter setting, p<n
c <- 0.4
n <- 1e2
p <- c*n
gamma <- 1

#### Parameter setting, p>n
c_g <- 1.2
n_g <- 1e2
p_g <- c_g*n_g

#### Generate samples of asset returns
set.seed(2)
X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)

set.seed(2)
X_g <- matrix(data=rnorm(n_g*p_g), nrow=p_g, ncol=n_g)


#### tests ####

# Interface -------------------------------------------------------------------

test_that("The interface of MVShrinkPortfolio works correctly", {

  expect_error(MVShrinkPortfolio(x=X, gamma='text instead of numeric',
                                 type='shrinkage', b=b, beta = 0.05))
  expect_error(MVShrinkPortfolio(x=X, gamma=gamma,
                                 type=4, b=b, beta = 0.05))
  expect_error(MVShrinkPortfolio(x=3, gamma=gamma,
                                 type='shrinkage', b=b, beta = 0.05))
  expect_error(MVShrinkPortfolio(x=X, gamma=gamma,
                                 type='invalid type', b=b, beta = 0.05))
})

# Traditional portfolios ------------------------------------------------------

test_that("Traditional portfolios work properly, p<n", {

  trad_port <- MVShrinkPortfolio(x=X, gamma=gamma, type='traditional')

  expect_type(trad_port$call, "language")
  expect_true(inherits(trad_port$cov_mtrx, "matrix"))
  expect_true(inherits(trad_port$inv_cov_mtrx, "matrix"))
  expect_type(trad_port$means, "double")
  expect_type(trad_port$weights, "double")

  expect_type(trad_port$Port_Var, "double")
  expect_type(trad_port$Port_mean_return, "double")
  expect_type(trad_port$Sharpe, "double")

  expect_s3_class(object=trad_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


test_that("Traditional portfolios work properly, p>n", {

  trad_port <- MVShrinkPortfolio(x=X_g, gamma=gamma, type='traditional')

  expect_type(trad_port$call, "language")
  expect_true(inherits(trad_port$cov_mtrx, "matrix"))
  expect_true(inherits(trad_port$inv_cov_mtrx, "matrix"))
  expect_type(trad_port$means, "double")
  expect_type(trad_port$weights, "double")

  expect_type(trad_port$Port_Var, "double")
  expect_type(trad_port$Port_mean_return, "double")
  expect_type(trad_port$Sharpe, "double")

  expect_s3_class(object=trad_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


# Shrinkage portfolios --------------------------------------------------------

b <- rep(1/p,p)
b_g <- rep(1/p_g, p_g)

test_that("Shrinkage portfolios work properly, p<n", {

  shrink_port <- MVShrinkPortfolio(x=X, gamma=gamma, b=b,
                                   beta = 0.05, type='shrinkage')

  expect_type(shrink_port$call, "language")
  expect_true(inherits(shrink_port$cov_mtrx, "matrix"))
  expect_true(inherits(shrink_port$inv_cov_mtrx, "matrix"))
  expect_type(shrink_port$means, "double")
  expect_type(shrink_port$weights, "double")

  expect_type(shrink_port$Port_Var, "double")
  expect_type(shrink_port$Port_mean_return, "double")
  expect_type(shrink_port$Sharpe, "double")

  expect_s3_class(object=shrink_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


test_that("Shrinkage portfolios work properly, p>n", {

  shrink_port <- MVShrinkPortfolio(x=X_g, gamma=gamma, b=b_g,
                                   beta = 0.05, type='shrinkage')

  expect_type(shrink_port$call, "language")
  expect_true(inherits(shrink_port$cov_mtrx, "matrix"))
  expect_true(inherits(shrink_port$inv_cov_mtrx, "matrix"))
  expect_type(shrink_port$means, "double")
  expect_type(shrink_port$weights, "double")

  expect_type(shrink_port$Port_Var, "double")
  expect_type(shrink_port$Port_mean_return, "double")
  expect_type(shrink_port$Sharpe, "double")

  expect_s3_class(object=shrink_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


# Shrinkage portfolios, gamma=Inf ---------------------------------------------

test_that("Shrinkage portfolios work properly, p<n, gamma=Inf", {

  shrink_port <- MVShrinkPortfolio(x=X, gamma=Inf, b=b,
                                   beta = 0.05, type='shrinkage')

  expect_type(shrink_port$call, "language")
  expect_true(inherits(shrink_port$cov_mtrx, "matrix"))
  expect_true(inherits(shrink_port$inv_cov_mtrx, "matrix"))
  expect_type(shrink_port$means, "double")
  expect_type(shrink_port$weights, "double")

  expect_type(shrink_port$Port_Var, "double")
  expect_type(shrink_port$Port_mean_return, "double")
  expect_type(shrink_port$Sharpe, "double")

  expect_s3_class(object=shrink_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


test_that("Shrinkage portfolios work properly, p>n", {

  shrink_port <- MVShrinkPortfolio(x=X_g, gamma=Inf, b=b_g,
                                   beta = 0.05, type='shrinkage')

  expect_type(shrink_port$call, "language")
  expect_true(inherits(shrink_port$cov_mtrx, "matrix"))
  expect_true(inherits(shrink_port$inv_cov_mtrx, "matrix"))
  expect_type(shrink_port$means, "double")
  expect_type(shrink_port$weights, "double")

  expect_type(shrink_port$Port_Var, "double")
  expect_type(shrink_port$Port_mean_return, "double")
  expect_type(shrink_port$Sharpe, "double")

  expect_s3_class(object=shrink_port,
                  class="MeanVar_portfolio", exact = FALSE)
})


# The case p=n is not treated -------------------------------------------------

#### Parameter setting, p>n
n_e <- 1e2
p_e <- n_e
b <- rep(1/p_e, p_e)

#### Generate samples of asset returns
set.seed(2)
X_e <- matrix(data=rnorm(n_e*p_e), nrow=p_e, ncol=n_e)

test_that("The case p=n is not treated", {

  expect_error(MVShrinkPortfolio(x=X_e, gamma=gamma, b=b_e,
                                 beta = 0.05, type='shrinkage'))
  expect_error(MVShrinkPortfolio(x=X_e, gamma=Inf, b=b_e,
                                 beta = 0.05, type='shrinkage'))
})


