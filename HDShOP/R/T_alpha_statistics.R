
d_0 <- function(gamma, p, n){

  c((1+1/(1-p/n))/gamma,
    -1,
    1/(gamma^2*(1-p/n)),
    (-1-1/(1-p/n))/gamma,
    1)
}


#
# t(d_0) * Omega_hat_al_c * d_0

# This function contains both c and c_n. Is that right?
Omega_hat_al_c <- function(x, b){

  c_n <- nrow(x)/ncol(x)
  M <- matrix(data=rep(0,25), nrow=5, ncol=5)
  s_hat_c <- s_hat_c(x)
  V_hat_c <- V_hat_GMV(x)/(1-c_n)
  V_hat_b <- V_hat_b(x, b)

  diag(M) <- c(V_hat_c*(s_hat_c+1)/(1-c_n), (2*V_hat_c^2)/(1-c_n),
               2*((s_hat_c+1)^2+c_n-1)/(1-c_n), V_hat_b, 2*V_hat_b^2)

  R_hat_b <- R_hat_b(x=x, b=b)
  R_hat_GMV<-R_hat_GMV(x=x)

  M[4,1] <- M[1,4] <- V_hat_c
  M[5,1] <- M[1,5] <- -2*V_hat_c*(R_hat_b-R_hat_GMV)
  M[5,2] <- M[2,5] <- 2*V_hat_c^2
  M[4,3] <- M[3,4] <- 2*(R_hat_b-R_hat_GMV)
  M[5,3] <- M[3,5] <- -2*(R_hat_b-R_hat_GMV)^2
  M
}


#  T_alpha, formula (44) in BDOP20
#
#' Test for MV portfolio weights
#'
#' An asymptotic test of a given mean-variance portfolio in a high-dimensional setting.
#' The tested hypotheses are
#' \deqn{H_0: w_{MV} = w_0 \quad vs \quad H_1: w_{MV} \neq w_0.}
#' The test statistic is based on the shrinkage estimator of MV portfolio weights
#' \insertCite{@see Eq.44 of @BDOPS2021}{HDShOP}.
#'
#' Note: when gamma == Inf, we get the test for the weights of the global minimum
#' variance portfolio as in Theorem 2 of \insertCite{BDPS2019;textual}{HDShOP}.
#' @inheritParams MVShrinkPortfolio
#' @param w_0 a numeric vector of tested weights.
#' @param beta a confidence level for the test.
#' @return
#'
#' | Element | Description |
#' | --- | --- |
#' | alpha_hat | the estimated shrinkage intensity |
#' | alpha_sd | the standard deviation of the shrinkage intensity |
#' | alpha_lower | the lower bound for the shrinkage intensity |
#' | alpha_upper | the upper bound for the shrinkage intensity  |
#' | T_alpha | the value of the test statistic |
#' | p_value | the p-value for the test |
#' @md
#' @references \insertAllCited{}
#' @examples
#' n<-3e2 # number of realizations
#' p<-.5*n # number of assets
#' b<-rep(1/p,p)
#' gamma<-1
#'
#' x <- matrix(data = rnorm(n*p), nrow = p, ncol = n)
#'
#' T_alpha <- test_MVSP(gamma=gamma, x=x, w_0=b, beta=0.05)
#' T_alpha
#' @export
test_MVSP <- function(gamma, x, w_0, beta=0.05) {

  n <- ncol(x)
  p <- nrow(x)

  Omega_hat_al_c <- Omega_hat_al_c(x=x, b=w_0)
  d_0<-d_0(gamma, p, n)
  B_hat <- B_hat(gamma=gamma, x=x, b=w_0)

  alpha_hat<-alpha_hat_star_c(gamma=gamma, x=x, b=w_0)
  alpha_sd<-as.numeric(sqrt(t(d_0) %*% Omega_hat_al_c %*% d_0) / B_hat/sqrt(n))
  z<-qnorm(p=1-beta/2 , mean = 0, sd = 1)
  alpha_lower<-alpha_hat-z*alpha_sd
  alpha_upper<-alpha_hat+z*alpha_sd
  T_alpha <- as.numeric(alpha_hat/ alpha_sd)
  p_value <- 2*(1-pnorm(abs(T_alpha), mean = 0, sd = 1))

  list(alpha_hat=alpha_hat, alpha_sd=alpha_sd, alpha_lower=alpha_lower,
       alpha_upper=alpha_upper, T_alpha=T_alpha, p_value=p_value)
}


