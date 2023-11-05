context("Shrinkage estimators for covariance")

#### Parameter setting
c <- 0.7
n <- 5e2
p <- c*n
gamma <- 1

mu_n <- rep(0, p)
Sigma <- matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- solve(Sigma)

#### Generate a sample of asset returns
X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)

#### Compute shrinkage estimators
Cov_LW02 <- nonlin_shrinkLW(X)

SCM <- Sigma_sample_estimator(X)
Cov_BGP14 <- CovShrinkBGP14(n=n, TM=Sigma, SCM=SCM)$S

#### tests
test_that("All cov estimators produce matrices", {

  expect_is(Cov_LW02, class="matrix")
  expect_is(Cov_BGP14, class="matrix")
})

test_that("All cov estimators are close to the real covariance", {

  expect_lte(sum(abs(Cov_LW02-Sigma)[1:6,1:6]), 6/100*5)
  expect_lte(sum(abs(Cov_BGP14-Sigma)[1:6,1:6]), 6/100*5)
})

test_that("nonlin_shrinkLW works when c>1", {

  n<-100; c<-4.2; p<-c*n; mu <- rep(0, p)

  set.seed(1)
  Sigma <- HDShOP::RandCovMtrx(p=p)
  set.seed(1)
  X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))

  LW <- nonlin_shrinkLW(X)
  min_eig <- min(eigen(LW)$values)

  expect_gte(min_eig, 0.05)
  expect_is(LW , class="matrix")

  #### the same with another c and seed

  n<-100; c<-2; p<-c*n; mu <- rep(0, p)

  set.seed(2)
  Sigma <- HDShOP::RandCovMtrx(p=p)
  set.seed(2)
  X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))

  LW <- nonlin_shrinkLW(X)
  min_eig <- min(eigen(LW)$values)

  expect_gte(min_eig, 0.05)
  expect_is(LW , class="matrix")
})

