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

b <- c(0.1,0.2,0.2,0.2,0.3)
test_that("R's are equivalent in the simplest case", {
  expect_lt(abs(R_hat_GMV(x)-R_GMV(mu=rep(0,5), Sigma=Sigma)), 0.0007)
  expect_lt(abs(R_hat_b(x=x, b=b)), 0.0007)
})


n<-1e4 # number of realizations
mu <- c(-1,5,3,-2,8)
x <-matrix(data = rnorm(n*p, mean=mu), nrow = p, ncol = n)

test_that("R's are equivalent to their asymptotics", {
  expect_lt(abs(R_hat_b(x=x, b=b) - R_b(mu=mu, b=b)), 0.01)
  expect_lt(abs(R_GMV(mu, Sigma) - R_hat_GMV(x)), 0.2)
})

test_that("V's are equivalent to their asymptotics", {
  expect_lt(abs(V_hat_b(x=x, b=b) - V_b(Sigma=Sigma, b=b)), 0.01)
  expect_lt(abs(V_hat_GMV(x) - V_GMV(Sigma)), 0.05)
})

test_that("s's are equivalent to their asymptotics", {
  expect_lt(abs(s(mu, Sigma) - s_hat_c(x)), 25)
})


n <- 1e3
gamma<-1
c <- 0.5
p <- c*n
x <- rnorm(n=n*p, mean=0)
x <- matrix(data=x, nrow=p, ncol=n)
mu<- rep(0, p)
b <- rep(1/p, p)

Sigma <- matrix(data=0, nrow=p, ncol=p)
diag(Sigma) <- 1

test_that("alpha_hat is equivalent to its asymptotic value", {
  expect_lt(abs(alpha_star(gamma=gamma, mu=mu, Sigma=Sigma, b=b, c=c) - alpha_hat_star_c(gamma, x, b)), 0.2)
})









