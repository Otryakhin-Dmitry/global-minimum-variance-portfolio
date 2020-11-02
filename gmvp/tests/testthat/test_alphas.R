context("Different expressions for alphas")

#### alphas EU and GMV work and are equal when gamma=Inf

test_that("alphas EU and GMV work and are equal when gamma=Inf", {

  n<-1e3 # number of realizations
  p<-0.8*n # number of assets
  w_0 <- rep(1/p,p)

  x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)

  al_EU <- alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
  al_GMV <- alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous

  expect_is(al_EU, class="numeric")
  expect_is(al_GMV, class="numeric")
  expect_equal(al_EU, al_GMV)
})


#### x with dependent components

test_that("alphas EU and GMV work and are equal when gamma=Inf and components of x are dependent", {

  if (!requireNamespace("MASS", quietly =TRUE)) skip("package MASS is not installed")
  library('MASS')

  n<-3e2 # number of realizations
  p<-0.3*n # number of assets
  w_0 <- rep(1/p,p)
  mu <- seq(0.2,-0.2, length.out=p)

  cov_mat <- SRandCovMtrx(n=n, p=p, mu=mu)
  x <- t(mvrnorm(n,mu, cov_mat))

  al_EU <- alpha_hat_star_c(gamma=Inf, x, b=w_0) # must be computable
  al_GMV <- alpha_hat_star_c_GMV(x, b=w_0) # must be equal to the previous

  expect_is(al_EU, class="numeric")
  expect_is(al_GMV, class="numeric")
  expect_equal(al_EU, al_GMV)
})


####
test_that("Deterministic alphas EU and GMV work and are equal when gamma=Inf", {

  n<-3e2 # number of realizations
  p<-0.3*n # number of assets
  w_0 <- rep(1/p,p)

  Sigma <- matrix(0, p, p)
  diag(Sigma) <- 1
  a_s_EU <- alpha_star(gamma=Inf, mu=rep(0,p), Sigma=Sigma, c=p/n, b=w_0) # must be computable
  a_s_GMV <- alpha_star_GMV(Sigma=Sigma, c=p/n, b=w_0) # must be equal to the previous

  expect_is(a_s_EU, class="numeric")
  expect_is(a_s_GMV, class="numeric")
  expect_equal(a_s_EU, a_s_GMV)

})





