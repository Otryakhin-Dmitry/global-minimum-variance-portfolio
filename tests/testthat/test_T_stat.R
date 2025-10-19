# context("T_alpha test")

###############################
n <- 5e2 # number of realizations
p <- 0.3*n # number of assets
w_0 <- rep(1/p,p)
gamma <- 5e-3
###############################


test_that("Ts from (44) and (41) are equivalent", {

  if (!requireNamespace("waldo", quietly =TRUE)) skip("package waldo is not installed")
  library('waldo')

  mu <- rep(0,p)
  d_0 <- d_0(gamma=gamma, p=p, n=n)
  ind_len <- 5e1
  X <- matrix(data = NA, nrow = ind_len, ncol = 2)

  Sigma <- matrix(data=0, nrow=p, ncol=p)
  diag(Sigma) <- 1

  for (ind in 1:ind_len) {

    x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
    X[ind,1] <- test_MVSP(gamma=gamma, x=x, w_0=w_0)$T_alpha

    # Computing T_alpha in another way, (29 & 41)
    Omega.est <- Omega_hat_al_c(x=x, b=w_0)

    t.a <- c(R_hat_GMV(x)-R_GMV(mu, Sigma),
            V_hat_c(x)-V_GMV(Sigma),
            s_hat_c(x)-s(mu, Sigma),
            R_hat_b(x, b=w_0)-R_b(mu, b=w_0),
            V_hat_b(x, b=w_0)-V_b(Sigma, b=w_0))

    X[ind,2]<- sqrt(n)*t(d_0)%*%t.a/sqrt(t(d_0)%*% Omega.est%*%d_0)
  }

  expect_true(all(round(X[,1], digits =4)==round(X[,2], digits =4)))
})



test_that("Ts are standard normal", {

  if (!requireNamespace("EstimDiagnostics", quietly =TRUE)) skip("package EstimDiagnostics is not installed")
  library('EstimDiagnostics')

  set.seed(1)
  sample <-
  replicate(n=3e1, {
    x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
    T_alpha <- test_MVSP(gamma=gamma, x=x, w_0=w_0, beta=0.05)
    T_alpha$T_alpha
  })
  expect_distfit(sample, p_value = 0.01, nulldist=pnorm, mean = 0, sd = 1)

  set.seed(2)
  sample <-
    replicate(n=3e1, {
      x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
      T_alpha <- test_MVSP(gamma=gamma, x=x, w_0=w_0, beta=0.05)
      T_alpha$T_alpha
    })
  expect_distfit(sample, p_value = 0.01, nulldist=pnorm, mean = 0, sd = 1)
})


