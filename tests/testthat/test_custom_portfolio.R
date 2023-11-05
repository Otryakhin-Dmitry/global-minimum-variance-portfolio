
#---------------------------------------------------------------

n<-100 
c<-.8 
p<-c*n 
mu <- rep(0, p)

set.seed(1)
Sigma <- HDShOP::RandCovMtrx(p=p)
set.seed(1)
X <- t(MASS::mvrnorm(n=n, mu=mu, Sigma=Sigma))

cov_mtrx <- CovarEstim(X, type="LW20", TM=TM)
means <- rowMeans(X)


#---------------------------------------------------------------

test_that("new_MeanVar_portfolio works properly", {

  cust_port_BS_LW <- new_MeanVar_portfolio(mean_vec=means,
                                           cov_mtrx=cov_mtrx, gamma=2)

  expect_type(cust_port_BS_LW$call, "language")
  expect_is(cust_port_BS_LW$cov_mtrx, "matrix")
  expect_is(cust_port_BS_LW$inv_cov_mtrx, "matrix")
  expect_type(cust_port_BS_LW$means, "double")
  expect_type(cust_port_BS_LW$weights, "double")

  expect_type(cust_port_BS_LW$Port_Var, "double")
  expect_type(cust_port_BS_LW$Port_mean_return, "double")
  expect_type(cust_port_BS_LW$Sharpe, "double")
  expect_s3_class(object=cust_port_BS_LW,
                  class="MeanVar_portfolio")

  expect_no_error(validate_MeanVar_portfolio(cust_port_BS_LW))
})


#---------------------------------------------------------------

test_that("MeanVar_portfolio works properly", {

  cust_port_BS_LW <- MeanVar_portfolio(mean_vec=means,
                                       cov_mtrx=cov_mtrx, gamma=2)

  expect_type(cust_port_BS_LW$call, "language")
  expect_is(cust_port_BS_LW$cov_mtrx, "matrix")
  expect_is(cust_port_BS_LW$inv_cov_mtrx, "matrix")
  expect_type(cust_port_BS_LW$means, "double")
  expect_type(cust_port_BS_LW$weights, "double")

  expect_type(cust_port_BS_LW$Port_Var, "double")
  expect_type(cust_port_BS_LW$Port_mean_return, "double")
  expect_type(cust_port_BS_LW$Sharpe, "double")
  expect_s3_class(object=cust_port_BS_LW,
                  class="MeanVar_portfolio")

  expect_no_error(validate_MeanVar_portfolio(cust_port_BS_LW))
})


#---------------------------------------------------------------

test_that("portfolio methods works properly", {

  cust_port_BS_LW <- MeanVar_portfolio(mean_vec=means,
                                       cov_mtrx=cov_mtrx, gamma=2)
  expect_no_error(summary(cust_port_BS_LW))
  expect_no_error(plot(cust_port_BS_LW))
})


#---------------------------------------------------------------

test_that("portfolio validator works properly", {

  cust_port_BS_LW <- new_MeanVar_portfolio(mean_vec=means,
                                           cov_mtrx=cov_mtrx, gamma=2)
  cust_port_BS_LW_no_cov <- cust_port_BS_LW
  cust_port_BS_LW_no_cov$cov_mtrx <- NA
  expect_error(validate_MeanVar_portfolio(cust_port_BS_LW_no_cov))

  cust_port_BS_LW_no_inv_cov <- cust_port_BS_LW
  cust_port_BS_LW_no_inv_cov$inv_cov_mtrx <- NA
  expect_error(validate_MeanVar_portfolio(cust_port_BS_LW_no_inv_cov))

  cust_port_BS_LW_no_means <- cust_port_BS_LW
  cust_port_BS_LW_no_means$means <- NA
  expect_error(validate_MeanVar_portfolio(cust_port_BS_LW_no_means))

  cust_port_BS_LW_no_weights <- cust_port_BS_LW
  cust_port_BS_LW_no_weights$weights <- NA
  expect_error(validate_MeanVar_portfolio(cust_port_BS_LW_no_weights))

})


