context("Work of various quantities")

###############################
p<-5 # number of assets
n<-1e7 # number of realizations
###############################

sm <- c(1,0,0,0,0,
        0,1,0,0,0,
        0,0,1,0,0,
        0,0,0,1,0,
        0,0,0,0,1)
Sigma <- matrix(data = sm, nrow = p, ncol = p)
x <-matrix(data = rnorm(n*p), nrow = p, ncol = n)

test_that("Non-diagonal elements are small when components are independent", {
  expect_lt(abs(R_hat_GMV(x)-R_GMV(mu=rep(0,5), Sigma=Sigma)), 0.0007)
})


b <- c(0.1,0.2,0.2,0.2,0.3)
test_that("R's are equivalent", {
  expect_lt(abs(R_hat_b(x=x, b=b)), 0.0007)
})


n<-1e3 # number of realizations
mu <- c(-1,5,3,-2,8)

x <-matrix(data = rnorm(n*p, mean=mu), nrow = p, ncol = n)

test_that("R's and V's are equivalent", {
  expect_lt(abs(R_hat_b(x=x, b=b) - R_b(mu=mu, b=b)), 0.6)
  expect_lt(abs(V_hat_b(x=x, b=b) - V_b(Sigma=Sigma, b=b)), 0.03)
})




