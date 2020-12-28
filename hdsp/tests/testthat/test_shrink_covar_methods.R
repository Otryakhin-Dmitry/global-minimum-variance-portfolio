context("Shrinkage estimators for covariance")

c <- 0.5
n <- 1e3
p <- c*n
gamma <- 1

mu_n<-rep(0, p)
Sigma=matrix(0, p, p)
diag(Sigma) <- 1
invSigma <- solve(Sigma)

X <- matrix(data=rnorm(n*p), nrow=p, ncol=n)

y_n_aver <- rowMeans(X)
Sigma_n_inv <- solve(Sigma)



Cov_LW02 <- nonlin_shrinkLW(X)



test_that("All cov estimators produce matrices", {

  expect_is(Cov_LW02, class="numeric")

})


test_that("All cov estimators are close to the real covariance", {

  expect_lte(sum(abs(Cov_LW02-Sigma)[1:6,1:6]), 36/100*5)

})











