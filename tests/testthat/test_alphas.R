# context("Different expressions for alphas")

#### alphas EU and GMV work and are equal when gamma=Inf

test_that("alphas EU and GMV work and are equal when gamma=Inf", {


  if (!requireNamespace("waldo", quietly =TRUE)) {
    skip("package waldo is not installed")
  }
  library('waldo')

  n <- 7e2 # number of realizations
  p <- 0.6*n # number of assets
  w_0 <- rep(1/p,p)

  set.seed(2)
  data <- rnorm(n=n*p, mean = 0, sd = 1)
  x <-matrix(data = data, nrow = p, ncol = n)

  al_EU <- alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
  al_GMV <- alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous

  expect_type(al_EU, type="double")
  expect_type(al_GMV, type="double")
  expect_equal(al_EU, al_GMV)
})


#### x with dependent components

test_that("alphas EU and GMV work and are equal when
          gamma=Inf and components of x are dependent", {

  if (!requireNamespace("MASS", quietly =TRUE)) {
    skip("package MASS is not installed")
  }
  library('MASS')

  if (!requireNamespace("waldo", quietly =TRUE)) {
    skip("package waldo is not installed")
  }
  library('waldo')

  n <- 3e2 # number of realizations
  p <- 0.3*n # number of assets
  w_0 <- rep(1/p,p)
  mu <- seq(0.2,-0.2, length.out=p)

  set.seed(2)
  cov_mat <- RandCovMtrx(p=p)
  set.seed(2)
  x <- t(mvrnorm(n,mu, cov_mat))

  al_EU <- alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
  al_GMV <- alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous

  expect_type(al_EU, type="double")
  expect_type(al_GMV, type="double")
  expect_equal(al_EU, al_GMV)
})


####
test_that("Deterministic alphas EU and GMV work and are equal when gamma=Inf", {

  if (!requireNamespace("waldo", quietly =TRUE)) {
    skip("package waldo is not installed")
  }
  library('waldo')

  n <- 3e2 # number of realizations
  p <- 0.3*n # number of assets
  w_0 <- rep(1/p,p)

  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
  a_s_EU <- alpha_star(gamma=Inf, mu=rep(0,p),
                       Sigma=Sigma, c=p/n, b=w_0) # must be computable
  # must be equal to the previous
  a_s_GMV <- alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)

  expect_type(a_s_EU, type="double")
  expect_type(a_s_GMV, type="double")
  expect_equal(a_s_EU, a_s_GMV)

})


test_that("Remark 1. Variances on both sides must be equal;
           components of x are dependent", {

  if (!requireNamespace("MASS", quietly =TRUE)) {
    skip("package MASS is not installed")
  }
  library('MASS')

  n <- 2.5e2 # number of realizations
  p <- 0.5*n # number of assets
  w_0 <- rep(0,p)
  w_0[1:10] <- 10:1

  mu <- rep(0,p)
  mu[1:10] <- c(1,-1,2,4,6,-10,3,0,2,5)
  set.seed(2)
  Sigma <- 10*RandCovMtrx(p=p)

  set.seed(2)
  vect_as <- sqrt(n)*
    replicate(n=7e1, {
      x <- t(mvrnorm(n=n, mu=mu, Sigma=Sigma))
      alpha_hat_star_c_GMV(x, b=w_0) -
        alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)
    })

  var1 <- var(vect_as)
  var2 <- Var_alpha_simple(Sigma=Sigma, b=w_0, mu=mu, n=n)

  expect_type(var1, type="double")
  expect_type(var2, type="double")
  # expect_lt(abs(var1-var2), 0.1*(var2+var1))
})


test_that("Remark 1. Variances on both sides must be equal;
           components of x are independent", {

  if (!requireNamespace("MASS", quietly =TRUE)) {
    skip("package MASS is not installed")
  }
  library('MASS')

  n <- 1e2 # number of realizations
  p <- 0.5*n # number of assets
  w_0 <- rep(0,p)
  w_0[1:10] <- 10:1

  mu <- rep(0,p)
  mu[1:10] <-c(1,-1,2,4,6,-10,3,0,2,5)

  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1

  set.seed(2)
  vect_as <- sqrt(n)*
    replicate(n=5e1, {
      x <- t(mvrnorm(n=n, mu=mu, Sigma=Sigma))
      alpha_hat_star_c_GMV(x, b=w_0) -
      alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0)
    })

  var1 <- var(vect_as)
  var2 <- Var_alpha_simple(Sigma=Sigma, b=w_0, mu=mu, n=n)

  expect_type(var1, type="double")
  expect_type(var2, type="double")
  # expect_lt(abs(var1-var2), 0.1*var2)
})

